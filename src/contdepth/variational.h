#ifndef VARIATIONAL_H
#define VARIATIONAL_H

#include "autorecon/stereo/flow_variational.h"
#include "autorecon/stereo/multigrid.h"
#include "EpipolarLine.h"

#include <nmisc/timer.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>

#include <vector>


/** 
    Warping and then taking the gradient is not entirely correct!!!!

    If we want more images, then the data term is sum over the phi's.
    
    Color channels.
 */
template <bool UseDisparity>
class VariationalDispProblem {
public:
  typedef EpipolarLineTemplate<UseDisparity> EpipolarLineType;

  constexpr static double omega = 1.25;
  mg::restriction_image restrict;

  double beta;  //!<Weight of smoothness term.
  std::vector<double>   S_alpha; // !<Weight of data terms
  std::vector<nacb::Imagef> S_I; // !<Data term tensors

  nacb::Imagef D_0, d;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  nacb::Imagef weight; //!<Per-pixel weighting.

  double              hx, hy;
  
  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);


  VariationalDispProblem(int width, int height){ 
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;    
  }

  VariationalDispProblem(const nacb::Imagef & im1, 
			 const std::vector<EpipolarLineTemplate<UseDisparity> > & elines, 
			 const std::vector<nacb::Imagef> & im2s,
			 const nacb::Imagef & _D_0,
			 double _alpha = 1.0,
			 double _beta = 6.0/255.0,
			 double _hx = 1.0, double _hy = 1.0) : 
    beta(_beta), D_0(_D_0), hx(_hx), hy(_hy) {

    assert(elines.size() == im2s.size());
    
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;

    f_h = nacb::Imagef(im1.w, im1.h, 1);
    d = nacb::Imagef(D_0.w, D_0.h, 1);
    residuals = nacb::Imagef(im1.w, im1.h, 1);
    weight = nacb::Imagef(im1.w, im1.h, 1);
    weight = 1;

    // Set the weights using the input image.
    if(im1.nchannels == 4){
      for(int y=0; y<im1.h; y++){
	for(int x=0; x<im1.w; x++){
	  if(im1(x, y, 3) < 0.1)
	    weight(x, y) = 0;
	}
      }
      weight = weight.boxFilter(2);
    }

    for(int i=0; i<(int)elines.size(); i++){
      S_I.push_back(nacb::Imagef(im1.w, im1.h, 3));
      S_alpha.push_back(_alpha);
      
      S_I[i] = 0;
    }

    // A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    d = 0;
    f_h = 0;

    int width = im1.w, height = im1.h;

    for(int i=0; i<(int)elines.size(); i++){
      const EpipolarLineType & eline = elines[i];
      const nacb::Imagef & im2 = im2s[i];
      const int nchans = std::min(3, im2.nchannels);
      
      // Assume that im2 is not warped
      nacb::Imagef im2_grad[3] = {im2.gradient(0),
				  im2.gradient(1),
				  im2.gradient(2)};

      nacb::Imagef Ix[3] = {im2_grad[0].getChannel(0)*(1.0/hx),
			    im2_grad[1].getChannel(0)*(1.0/hx),
			    im2_grad[2].getChannel(0)*(1.0/hx)};
			    
      nacb::Imagef Iy[3] = {im2_grad[0].getChannel(1)*(1.0/hy),
			    im2_grad[1].getChannel(1)*(1.0/hy),
			    im2_grad[2].getChannel(1)*(1.0/hy)};

      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  nacb::Vec2f co = eline(x, y, _D_0(x, y));
	  nacb::Vec2f co_deriv = eline.disparityDeriv(x, y, _D_0(x, y));
	
	  for(int chan = 0; chan<nchans; chan++){
	    double der = Ix[chan].bilinear(co.x, co.y)*co_deriv.x + 
	      Iy[chan].bilinear(co.x, co.y)*co_deriv.y;
	    
	    double Iz = im2.bilinear(co.x, co.y, chan) - im1(x, y, chan);
	    
	    // Image data term
	    S_I[i](x, y, 0) += der*der;
	    S_I[i](x, y, 1) += der*Iz;
	    S_I[i](x, y, 2) += Iz*Iz;
	  }
	}
      }
    }
  }
  

  nacb::Vec2f forwardBackward(const EpipolarLineType & eline,
			      const EpipolarLineType & eline_inv,
			      const nacb::Imagef & disp,
			      float x, float y, float d){
    nacb::Vec2f co = eline(x, y, d);
    float dback = disp.bilinear(co.x, co.y);
    return eline_inv(co.x, co.y, dback) - nacb::Vec2f(x, y);    
  }

    //Doesn't work ....
  void addDisplacementConsistency(const std::vector<EpipolarLineType> & elines,
				  const std::vector<EpipolarLineType> & eline_invs,
				  const std::vector<nacb::Imagef> & disp2s,
				  double _alpha = 0.01){
    int S_base = S_I.size();

    // Add the extra terms.
    int width = f_h.w, height = f_h.h;
    for(int i=0; i<(int)elines.size(); i++){
      S_I.push_back(nacb::Imagef(width, height, 3));
      S_alpha.push_back(_alpha);
      
      S_I.back() = 0;
    }

    nacb::Imagef res(width, height, 1);
    res = 0;

    double totalDiff = 0;
    int ndiff = 0;

    // Compute the actual derivatives.
    for(int i=0; i<(int)elines.size(); i++){
      const EpipolarLineType & eline = elines[i];
      const nacb::Imagef & disp2 = disp2s[i];
      
      // Assume that im2 is not warped
      nacb::Imagef d2_grad = disp2.gradient(0);
      nacb::Imagef Ix = d2_grad.getChannel(0)*(1.0/hx);
      nacb::Imagef Iy = d2_grad.getChannel(1)*(1.0/hy);

      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  if(weight(x, y) >= 1.0){
	    nacb::Vec2f co = eline(x, y, D_0(x, y));
	    nacb::Vec2f co_deriv = eline.disparityDeriv(x, y, D_0(x, y));
	    	    
	    float codisp = disp2.bilinear(co.x, co.y);
	    nacb::Vec2f co2 = eline_invs[i](co.x, co.y, codisp);
	    nacb::Vec2f co2_deriv = eline_invs[i].disparityDeriv(co.x, co.y, codisp);

	    nacb::Vec2f co2_deriv_xy[2];
	    eline_invs[i].deriv_xy(co.x, co.y, codisp, co2_deriv_xy);
	    
	    // Partial co2 partial disp
	    nacb::Vec2f og(x, y);

#ifdef  CHECK_FD
	    double eps = 1e-4;
	    nacb::Vec2f f = forwardBackward(elines[i], eline_invs[i], disp2, x, y, D_0(x, y));
	    nacb::Vec2f fup = forwardBackward(elines[i], eline_invs[i], disp2, x, y, D_0(x, y) + eps);
	    nacb::Vec2f fdn = forwardBackward(elines[i], eline_invs[i], disp2, x, y, D_0(x, y) - eps);

	    fup = (fup - fdn)*(0.5/eps);
	    
	    //printf("-----\n");
	    //printf("%f %f\n", f.x, f.y);
	    //printf("%f %f\n", fup.x, fup.y);
#endif


	    for(int k = 0; k<2; k++){
	      double der = co2_deriv_xy[k].data[0]*co_deriv.x +
		co2_deriv_xy[k].data[1]*co_deriv.y + 
		co2_deriv.data[k]*(Ix.bilinear(co.x, co.y)*co_deriv.x +
				   Iy.bilinear(co.x, co.y)*co_deriv.y);
	      
	      double Iz = co2.data[k] - og.data[k];
	      
#ifdef  CHECK_FD
	      //printf("%f %f %f\n", Iz, der, fup.data[k]);
	      //der = fup.data[k];
#endif
	      
	      res(x, y) += Iz*Iz;
	      totalDiff += Iz*Iz;
	      ndiff++;
	     
	      /*
	      if(fabs(der)>10000)
		printf("%f\n", der);
	      if(der>100)
		der = 100;
	      if(der<-100)
		der = -100;
	      */
	      if(isnan(Iz) || isnan(der) || isinf(Iz) || isinf(der))
		continue;

	      if(fabs(Iz) > 1e4 || fabs(der) > 1e4)
		continue;

	      // Image data term
	      S_I[i + S_base](x, y, 0) += der*der;
	      S_I[i + S_base](x, y, 1) += der*Iz;
	      S_I[i + S_base](x, y, 2) += Iz*Iz;
	    }
	  }
	}
      }
    }    

    double n, x;
    res.getRange(n, x);
    printf("total res %f %f - %f\n", totalDiff/ndiff, n, x);
    static int i = 0;
    char fname[1024];
    snprintf(fname, 1024, "/tmp/r-%d.png", i);
    res.getNormalizedImage().write(fname);
    i++;

  }

  VariationalDispProblem  restriction();
  
  nacb::Imagef prolong(const nacb::Imagef & img){
    return restrict.prolong(img);
  }
  
  robust_stats iteration(nacb::Imagef & residuals, bool update = true);
   
  void smooth(){
    iteration(residuals, true);
  }

  void setDataWeight(const nacb::Imagef & weights){
    weight = weights.copy();
  }

  void add_to_solution(const nacb::Imagef & update){d += update;}
  nacb::Imagef & get_solution(){return d;}
  nacb::Imagef copy_solution(){return d.copy();}

  int get_width() const {return d.getWidth();}
  int get_height() const {return d.getHeight();}
  bool is_linear() const {return false;}
  bool can_restrict() const {return get_width()>=16 && get_height()>=16;}

  nacb::Imagef getCombinedResult(){
    return D_0 + d;
  }

  /**
     The following functions are for computing cost terms
     given a current solution.  Not useful in the context of the linearization,
     more useful for checking between successive iterations of the method.
  */
  static double imageDataEnergy(const nacb::Imagef & D,
				const nacb::Imagef & image,
				const nacb::Imagef & weight,
				const std::vector<EpipolarLineType> & elines, 
				const std::vector<nacb::Imagef> & im2s,
				double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	
	nacb::Vec3f color(image(x, y, 0), image(x, y, 1), image(x, y, 2));
	
	for(int si=0; si<(int)elines.size(); si++){
	  nacb::Vec2f co2 = elines[si](x, y, D(x, y));
	  nacb::Vec3f color2;
	  
	  im2s[si].bilinear(co2.x, co2.y, color2.data, 3);
	  
	  nacb::Vec3f diff = color - color2;
	  energy_D += weight(x, y)*phi_D(diff.dot(diff));
	}	
      }
    }
    return energy_D;
  }

  static double smoothnessEnergy(const nacb::Imagef & D,
				 double (* phi_S)(double ) = phi_1en6){
    double energy_S = 0;
    int width = D.width;
    int height = D.height;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	int xup = std::min(x + 1, width - 1);
	int yup = std::min(y + 1, height - 1);
	int xdn = std::max(x - 1, 0);
	int ydn = std::max(y - 1, 0);
	
	nacb::Vec2d grad((D(xup, y) - D(xdn, y))/2.0,
			 ((D(x, yup) - D(x, ydn))/2.0));
	
	energy_S += phi_S(grad.dot(grad));
      }
    }
    return energy_S;
  }

  static double geometryDataEnergy(const nacb::Imagef & D,
				   const nacb::Imagef & weight,
				   const std::vector<EpipolarLineType> & elines, 
				   const std::vector<EpipolarLineType> & eline_invs,
				   const std::vector<nacb::Imagef> & disps,
				   double (* phi_D)(double ) = phi_1en6){
    double energy = 0;
    int width = D.width;
    int height = D.height;

    assert(elines.size() == disps.size() && elines.size() == eline_invs.size());

    assert(D.w == weight.w && D.height == weight.h);
    
    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	
	for(int i=0; i<(int)elines.size(); i++){
	  nacb::Vec2f co = elines[i](x, y, D(x, y));

	  co.x = std::max(0.f, std::min(co.x, (float)width - 1));
	  co.y = std::max(0.f, std::min(co.y, (float)height - 1));

	  nacb::Vec2f co_back = eline_invs[i](co.x, co.y, disps[i].bilinear(co.x, co.y));
	  nacb::Vec2f diff = nacb::Vec2f(x, y) - co_back;
	  
	  energy += phi_D(diff.dot(diff))*weight(x, y);
	}
      }
    }
    return energy;
  }
};



#endif
