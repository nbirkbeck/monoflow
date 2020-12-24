#ifndef VARIATIONAL_H
#define VARIATIONAL_H

#include "DisplaceUV.h"
#include "vec5.h"
#include "Math.h"

#include "autorecon/stereo/flow_variational.h"
#include "autorecon/stereo/multigrid.h"

#include <nmisc/timer.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>

#include <nimage/image.h>

#include <vector>

using namespace nacb;

/**
 */
class VariationalOffsetProblemUV {
public:

  constexpr static double omega = 1.25;
  mg::restriction_image restrict;

  static Vec2d defaultBeta;
  static double defaultAlpha;

  Vec2d beta;  //!<Weight of smoothness term.
  std::vector<double>   S_alpha; // !<Weight of data terms
  std::vector<nacb::Imagef> S_I; // !<Data term tensors

  nacb::Imagef D_0, d;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  std::vector<nacb::Imagef> weights; //!<Per-pixel weighting, one for each S_I

  double              hx, hy;
  
  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);


  VariationalOffsetProblemUV(int width, int height){ 
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;    
  }

  VariationalOffsetProblemUV(const std::vector<std::vector<nacb::Imagef> > & images, 
			     const std::vector<DisplaceOffsetUV> & warps, // Mapping from UV-space to image
			     const std::vector<nacb::Imagef> & indyWeights,
			     std::vector<std::pair<int, int> > & pairs,
			     const nacb::Imagef & _D_0,
			     const Vec2d & _alpha = Vec2d(1, 1),
			     const Vec2d & _beta = VariationalOffsetProblemUV::defaultBeta,
			     double _hx = 1.0, double _hy = 1.0) : 
    beta(_beta), D_0(_D_0), hx(_hx), hy(_hy) {

    int width = D_0.w, height = D_0.h;
    
    assert(warps.size() == images.size());
    assert(indyWeights.size() == warps.size());
    
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;

    f_h = nacb::Imagef(width, height, 4);
    d = nacb::Imagef(width, height, 4);
    residuals = nacb::Imagef(width, height, 4);

    // The weighting for data terms is the multiplication of the two
    for (int i=0; i<(int)pairs.size(); i++) {
      weights.push_back(indyWeights[pairs[i].first] * 
			indyWeights[pairs[i].second]);
      weights.push_back(indyWeights[pairs[i].first] * 
			indyWeights[pairs[i].second]);
    }
    
    // The single weights.
    for (int i=0; i<(int)images.size(); i++) {
      weights.push_back(indyWeights[i]);
    }

    // Stereo terms, one for each time.
    for(int i=0; i<(int)pairs.size(); i++){
      for (int j=0; j<2; j++) {
	S_I.push_back(nacb::Imagef(width, height, 15));
	S_I.back() = 0;
	S_alpha.push_back(_alpha.x * 0.5); // Multiplied by 0.5 as there is one for t=0 and t=1
      }
    }

    int indyBase = S_I.size();
    for(int i=0; i<(int)images.size(); i++){
      S_I.push_back(nacb::Imagef(width, height, 15));
      S_I.back() = 0;
      S_alpha.push_back(_alpha.y);
    }

    // A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    d = 0;
    f_h = 0;
    

    for(int i=0; i<(int)pairs.size(); i++){
      printf("Linearizing around warps.\n");
      
      int index1 = pairs[i].first;
      int index2 = pairs[i].second;
      const DisplaceOffsetUV & warp1 = warps[index1];
      const DisplaceOffsetUV & warp2 = warps[index2];
      const std::vector<nacb::Imagef> & im1 = images[index1];
      const std::vector<nacb::Imagef> & im2 = images[index2];
      //{images[index1][0], images[index1][1]};
      //const nacb::Imagef & im2[2] = {images[index2][0], images[index2][1]};
      
      const int nchans = std::min(3, im2[0].nchannels);

      nacb::Imagef Ix1[2], Iy1[2];
      nacb::Imagef Ix2[2], Iy2[2];

      printf("Getting gradients %ld %ld.\n", im1.size(), im2.size());
      for (int k=0; k<2; k++) {
	printf("Gradient: %d  %dx%d  %dx%d\n", k, im1[k].w, im1[k].h, im2[k].w, im2[k].h);
	fflush(stdout);
	image_gradient_all(im1[k], Ix1[k], Iy1[k]);
	image_gradient_all(im2[k], Ix2[k], Iy2[k]);

	Ix1[k] *= (1.0/hx);
	Iy1[k] *= (1.0/hy);
	
	Ix2[k] *= (1.0/hx);
	Iy2[k] *= (1.0/hy);
      }

      printf("Evaluating.\n");
      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  Vec3f offs(D_0(x, y, 1), D_0(x, y, 2), D_0(x, y, 3));
	  nacb::Vec2d co1_o_deriv[3];
	  nacb::Vec2d co2_o_deriv[3];
	  

	  Vec3f zeroOffset(0, 0, 0);

	  nacb::Vec2f co1[2] = {warp1(x, y, D_0(x, y, 0), zeroOffset), warp1(x, y, D_0(x, y, 0), offs)};
	  nacb::Vec2f co2[2] = {warp2(x, y, D_0(x, y, 0), zeroOffset), warp2(x, y, D_0(x, y, 0), offs)};
	  nacb::Vec2d co1_deriv[3];
	  nacb::Vec2d co2_deriv[3];
	  
	  
	  co1_deriv[0] = warp1.derivatives(x, y, D_0(x, y, 0), zeroOffset);
	  co2_deriv[0] = warp2.derivatives(x, y, D_0(x, y, 0), zeroOffset);

	  co1_deriv[1] = warp1.derivatives(x, y, D_0(x, y, 0), offs, co1_o_deriv);
	  co2_deriv[1] = warp2.derivatives(x, y, D_0(x, y, 0), offs, co2_o_deriv);

	  for(int chan = 0; chan<nchans; chan++){
	    // Difference between two images at time t=0
	    {
	      Vec2d i1chan(Ix1[0].bilinear(co1[0].x, co1[0].y, chan),
			   Iy1[0].bilinear(co1[0].x, co1[0].y, chan));
	      
	      Vec2d i2chan(Ix2[0].bilinear(co2[0].x, co2[0].y, chan),
			   Iy2[0].bilinear(co2[0].x, co2[0].y, chan));
	      
	      double Iz = im2[0].bilinear(co2[0].x, co2[0].y, chan) - 
		im1[0].bilinear(co1[0].x, co1[0].y, chan);
	      
	      Vec5d II(i2chan.dot(co2_deriv[0]) - i1chan.dot(co1_deriv[0]),
		       0, 0, 0, // Not a function of offsets.
		       Iz);
	      
	      accumulate(S_I[2*i], x, y, Mat5x5::outerProduct(II));
	    }

	    Vec2d i1chan(Ix1[1].bilinear(co1[1].x, co1[1].y, chan),
			 Iy1[1].bilinear(co1[1].x, co1[1].y, chan));

	    Vec2d i2chan(Ix2[1].bilinear(co2[1].x, co2[1].y, chan),
			 Iy2[1].bilinear(co2[1].x, co2[1].y, chan));
	    
	    double Iz = im2[1].bilinear(co2[1].x, co2[1].y, chan) 
	      - im1[1].bilinear(co1[1].x, co1[1].y, chan);

	    Vec5d II(i2chan.dot(co2_deriv[1]) - i1chan.dot(co1_deriv[1]),
		     // These derivatives only exist for the t=1 case.
		     i2chan.dot(co2_o_deriv[0]) - i1chan.dot(co1_o_deriv[0]),
		     i2chan.dot(co2_o_deriv[1]) - i1chan.dot(co1_o_deriv[1]),
		     i2chan.dot(co2_o_deriv[2]) - i1chan.dot(co1_o_deriv[2]),
		     Iz);
	    // Difference between two images at time t=1
	    accumulate(S_I[2*i+1], x, y, Mat5x5::outerProduct(II));
	  }
	}
      }
    }

    for (int i=0; i<(int)images.size(); i++) {
      printf("Linearizing flow in time: %d\n", i);

      const DisplaceOffsetUV & warp = warps[i];
      const std::vector<nacb::Imagef> & im = images[i];
      const int nchans = std::min(3, im[0].nchannels);

      nacb::Imagef Ix[2], Iy[2];

      for (int k=0; k<2; k++) {
	image_gradient_all(im[k], Ix[k], Iy[k]);

	Ix[k] *= (1.0/hx);
	Iy[k] *= (1.0/hy);
      }

      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  Vec3f offs(D_0(x, y, 1), D_0(x, y, 2), D_0(x, y, 3));
	  nacb::Vec2d co_o_deriv[3];

	  Vec3f zeroOffset(0, 0, 0);

	  nacb::Vec2f co[2] = {warp(x, y, D_0(x, y, 0), zeroOffset), warp(x, y, D_0(x, y, 0), offs)};
	  nacb::Vec2d co_deriv[3];
	  co_deriv[0] = warp.derivatives(x, y, D_0(x, y, 0), zeroOffset);
	  co_deriv[1] = warp.derivatives(x, y, D_0(x, y, 0), offs, co_o_deriv);

	  for(int chan = 0; chan<nchans; chan++){
	    // Difference between images at I_1 and I_0
	    Vec2d t0chan(Ix[0].bilinear(co[0].x, co[0].y, chan),
			 Iy[0].bilinear(co[0].x, co[0].y, chan));
	      
	    Vec2d t1chan(Ix[1].bilinear(co[1].x, co[1].y, chan),
			 Iy[1].bilinear(co[1].x, co[1].y, chan));
	    
	    double Iz = im[1].bilinear(co[1].x, co[1].y, chan) - im[0].bilinear(co[0].x, co[0].y, chan);

	    Vec5d II(t1chan.dot(co_deriv[1]) - t0chan.dot(co_deriv[0]),
		     t1chan.dot(co_o_deriv[0]),
		     t1chan.dot(co_o_deriv[1]),
		     t1chan.dot(co_o_deriv[2]),
		     Iz);

	    accumulate(S_I[indyBase + i], x, y, Mat5x5::outerProduct(II));
	  }
	}
      }
    }

  }
 
  VariationalOffsetProblemUV  restriction();
  
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
				const std::vector<std::vector<nacb::Imagef> > & images,
				const std::vector<DisplaceOffsetUV> & warps,
				const std::vector<nacb::Imagef> & weight,
				const std::vector<std::pair<int, int> > & pairs,
				double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){	
	nacb::Vec3f offs(D(x, y, 1), D(x, y, 2), D(x, y, 3));

	for (int i=0; i<(int)pairs.size(); i++) {
	  int i1 = pairs[i].first, i2 = pairs[i].second;

	  for (int t=0; t<=1; t++) {
	    nacb::Vec2f co1 = warps[i1](x, y, D(x, y), offs*t);
	    nacb::Vec2f co2 = warps[i2](x, y, D(x, y), offs*t);
	    nacb::Vec3f color1, color2;
	    
	    images[i1][t].bilinear(co1.x, co1.y, color1.data, 3);
	    images[i2][t].bilinear(co2.x, co2.y, color2.data, 3);
	    
	    nacb::Vec3f diff = color2 - color1;
	    energy_D += weight[i1](x, y)*weight[i2](x, y)*phi_D(diff.dot(diff));
	  }
	}

	for (int i=0; i<(int)images.size(); i++) {
	  nacb::Vec2f co1 = warps[i](x, y, D(x, y), offs*0);
	  nacb::Vec2f co2 = warps[i](x, y, D(x, y), offs);
	  nacb::Vec3f color1, color2;
	
	  images[i][0].bilinear(co1.x, co1.y, color1.data, 3);
	  images[i][1].bilinear(co2.x, co2.y, color2.data, 3);
	  
	  nacb::Vec3f diff = color2 - color1;
	  energy_D += weight[i](x, y)*phi_D(diff.dot(diff));
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
	// FIXME: doesn't take into account the flow.
	energy_S += phi_S(grad.dot(grad));
      }
    }
    return energy_S;
  }
};



#endif
