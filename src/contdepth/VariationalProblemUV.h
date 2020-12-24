#ifndef VARIATIONAL_H
#define VARIATIONAL_H

#include "DisplaceUV.h"

#include "autorecon/stereo/flow_variational.h"
#include "autorecon/stereo/multigrid.h"

#include <nmisc/timer.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>

#include <nimage/image.h>

#include <vector>


/**
 */
class VariationalProblemUV {
public:

  constexpr static double omega = 1.25;
  mg::restriction_image restrict;

  double beta;  //!<Weight of smoothness term.
  std::vector<double>   S_alpha; // !<Weight of data terms
  std::vector<nacb::Imagef> S_I; // !<Data term tensors

  nacb::Imagef D_0, d;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  std::vector<nacb::Imagef> weights; //!<Per-pixel weighting, one for each S_I

  double              hx, hy;
  
  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);


  VariationalProblemUV(int width, int height){ 
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;    
  }

  VariationalProblemUV(const std::vector<nacb::Imagef> & images, 
		       const std::vector<DisplaceUV> & warps, // Mapping from UV-space to image
		       const std::vector<nacb::Imagef> & indyWeights,
		       std::vector<std::pair<int, int> > & pairs,
		       const nacb::Imagef & _D_0,
		       double _alpha = 1.0,
		       double _beta = 6.0/255.0,
		       double _hx = 1.0, double _hy = 1.0) : 
    beta(_beta), D_0(_D_0), hx(_hx), hy(_hy) {

    int width = D_0.w, height = D_0.h;
    
    assert(warps.size() == images.size());
    assert(indyWeights.size() == warps.size());
    
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;

    f_h = nacb::Imagef(width, height, 1);
    d = nacb::Imagef(width, height, 1);
    residuals = nacb::Imagef(width, height, 1);

    // The weighting for data terms is the multiplication of the two
    for (int i=0; i<(int)pairs.size(); i++) {
      weights.push_back(indyWeights[pairs[i].first] * 
			indyWeights[pairs[i].second]);
    }

    for(int i=0; i<(int)pairs.size(); i++){
      S_I.push_back(nacb::Imagef(width, height, 3));
      S_alpha.push_back(_alpha);
      
      S_I[i] = 0;
    }

    // A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    d = 0;
    f_h = 0;
    

    for(int i=0; i<(int)pairs.size(); i++){
      int index1 = pairs[i].first;
      int index2 = pairs[i].second;
      const DisplaceUV & warp1 = warps[index1];
      const DisplaceUV & warp2 = warps[index2];
      const nacb::Imagef & im1 = images[index1];
      const nacb::Imagef & im2 = images[index2];
      
      const int nchans = std::min(3, im2.nchannels);
      
      // Assume that the images are not warped.
      nacb::Imagef im1_grad[3] = {im1.gradient(0),
				  im1.gradient(1),
				  im1.gradient(2)};

      nacb::Imagef im2_grad[3] = {im2.gradient(0),
				  im2.gradient(1),
				  im2.gradient(2)};

      nacb::Imagef Ix1[3] = {im1_grad[0].getChannel(0)*(1.0/hx),
			     im1_grad[1].getChannel(0)*(1.0/hx),
			     im1_grad[2].getChannel(0)*(1.0/hx)};
			    
      nacb::Imagef Iy1[3] = {im1_grad[0].getChannel(1)*(1.0/hy),
			     im1_grad[1].getChannel(1)*(1.0/hy),
			     im1_grad[2].getChannel(1)*(1.0/hy)};

      nacb::Imagef Ix2[3] = {im2_grad[0].getChannel(0)*(1.0/hx),
			     im2_grad[1].getChannel(0)*(1.0/hx),
			     im2_grad[2].getChannel(0)*(1.0/hx)};
			    
      nacb::Imagef Iy2[3] = {im2_grad[0].getChannel(1)*(1.0/hy),
			     im2_grad[1].getChannel(1)*(1.0/hy),
			     im2_grad[2].getChannel(1)*(1.0/hy)};

      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  nacb::Vec2f co1 = warp1(x, y, _D_0(x, y));
	  nacb::Vec2f co2 = warp2(x, y, _D_0(x, y));
	  nacb::Vec2f co1_deriv = warp1.depthDeriv(x, y, _D_0(x, y));
	  nacb::Vec2f co2_deriv = warp2.depthDeriv(x, y, _D_0(x, y));

	  for(int chan = 0; chan<nchans; chan++){
	    double der1 = 
	      Ix1[chan].bilinear(co1.x, co1.y)*co1_deriv.x + 
	      Iy1[chan].bilinear(co1.x, co1.y)*co1_deriv.y;

	    double der2 = 
	      Ix2[chan].bilinear(co2.x, co2.y)*co2_deriv.x + 
	      Iy2[chan].bilinear(co2.x, co2.y)*co2_deriv.y;
	    
	    double Iz = im2.bilinear(co2.x, co2.y, chan) - im1.bilinear(co1.x, co1.y, chan);

	    double der = der2 - der1;

	    // Image data term
	    S_I[i](x, y, 0) += der*der;
	    S_I[i](x, y, 1) += der*Iz;
	    S_I[i](x, y, 2) += Iz*Iz;
	  }
	}
      }
    }
  }
 
  VariationalProblemUV  restriction();
  
  nacb::Imagef prolong(const nacb::Imagef & img){
    return restrict.prolong(img);
  }
  
  robust_stats iteration(nacb::Imagef & residuals, bool update = true);
   
  void smooth(){
    iteration(residuals, true);
  }

  void setDataWeights(const std::vector<nacb::Imagef> & weightsIn){
    weights.clear();
    
    for (int i=0; i<(int)weightsIn.size(); i++) {
      weights.push_back(weightsIn[i].copy());
    }
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
				const std::vector<nacb::Imagef> & images,
				const std::vector<DisplaceUV> & warps,
				const std::vector<nacb::Imagef> & weight,
				const std::vector<std::pair<int, int> > & pairs,
				double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){	
	for (int i=0; i<(int)pairs.size(); i++) {
	  int i1 = pairs[i].first, i2 = pairs[i].second;
	  
	  nacb::Vec2f co1 = warps[i1](x, y, D(x, y));
	  nacb::Vec2f co2 = warps[i2](x, y, D(x, y));
	  nacb::Vec3f color1, color2;
	  
	  images[i1].bilinear(co1.x, co1.y, color1.data, 3);
	  images[i2].bilinear(co2.x, co2.y, color2.data, 3);
	  
	  nacb::Vec3f diff = color2 - color1;
	  energy_D += weight[i1](x, y)*weight[i2](x, y)*phi_D(diff.dot(diff));
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
};



#endif
