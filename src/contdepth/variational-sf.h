#ifndef VARIATIONAL_SF_H
#define VARIATIONAL_SF_H

#include "autorecon/stereo/flow_variational.h"
#include "autorecon/stereo/multigrid.h"
#include "EpipolarLine.h"
#include <nmisc/timer.h>
#include <vector>
#include <nmath/vec2.h>

#include "vec5.h"
#include "Math.h"
#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"

#include "DisparityOffset.h"

using namespace nacb;

namespace std {
  inline ostream & operator<<(ostream & ostr, const nacb::Vec2d & o){
    ostr << "[" << o.x << "," << o.y << "]";
    return ostr;
  }
};



/** 
    Warping and then taking the gradient is not entirely correct!!!!

    If we want more images, then the data term is sum over the phi's.
    
    Color channels.
 */
template <bool UseDisparity = true>
class VariationalSceneFlowProblem {
public:

  static nacb::Vec2d defaultBeta;
  static double defaultAlpha;

  constexpr static double omega = 1.25;
  mg::restriction_image restrict;

  Vec2d beta;  //!<Weight of smoothness term.
  std::vector<double>   S_alpha; // !<Weight of data terms
  std::vector<nacb::Imagef> S_I; // !<Data term tensors

  nacb::Imagef D_0, d;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  nacb::Imagef weight; // !<Per-pixel weighting.

  double              hx, hy;

  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);


  VariationalSceneFlowProblem(int width, int height){ 
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;    
  }

  VariationalSceneFlowProblem(const nacb::Imagef & im1, 
			      const std::vector<DisparityOffsetBase<UseDisparity> > & warps, 
			      const std::vector<nacb::Imagef> & im2s,
			      const nacb::Imagef & _D_0,
			      double _alpha = VariationalSceneFlowProblem::defaultAlpha,
			      const Vec2d & _beta = VariationalSceneFlowProblem::defaultBeta, 
			      double _hx = 1.0, double _hy = 1.0) : 
    beta(_beta), D_0(_D_0), hx(_hx), hy(_hy) {

    assert(warps.size() == im2s.size());
    
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;

    f_h = nacb::Imagef(im1.w, im1.h, 4);
    d = nacb::Imagef(D_0.w, D_0.h, 4);
    residuals = nacb::Imagef(im1.w, im1.h, 4);

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

    for(int i=0; i<(int)warps.size(); i++){
      S_I.push_back(nacb::Imagef(im1.w, im1.h, 15));
      S_alpha.push_back(_alpha);
      
      S_I[i] = 0;
    }

    //A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    d = 0;
    f_h = 0;

    int width = im1.w, height = im1.h;
    for(int i=0; i<(int)warps.size(); i++){  
      const nacb::Imagef & im2 = im2s[i];
      const int nchans = std::min(3, im2.nchannels);
      const DisparityOffsetBase<UseDisparity> & warp = warps[i];

      // Fix the gradient computation (could be more efficient).

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
	  nacb::Vec3d offs(_D_0(x, y, 1),
			   _D_0(x, y, 2),
			   _D_0(x, y, 3));			   
	  nacb::Vec2d co;
	  nacb::Vec2d co_deriv[4];
	  
	  co = warp(x, y, _D_0(x, y, 0), offs, co_deriv);
	
	  for(int chan = 0; chan<nchans; chan++){
	    Vec2d ichan(Ix[chan].bilinear(co.x, co.y),
			Iy[chan].bilinear(co.x, co.y));
	    
	    double Iz = im2.bilinear(co.x, co.y, chan) - im1(x, y, chan);
	    Vec5d II(ichan.dot(co_deriv[0]),
		     ichan.dot(co_deriv[1]),
		     ichan.dot(co_deriv[2]),
		     ichan.dot(co_deriv[3]),
		     Iz);
	    accumulate(S_I[i], x, y, Mat5x5::outerProduct(II));
	  }
	}
      }
    }
  }  
  

  /**
     \brief This constraint only really penalizes the flow from being
            close to the surfaces in the other view (e.g., it doesn't
	    use the flow in the previos flow).
  */
  Vec2f forwardBackwardSimple(const DisparityOffsetBase<UseDisparity> & warp,
			      const DisparityOffsetBase<UseDisparity> & warp_inv,
			      const Imagef & flow,
			      int x, int y,
			      double disp,
			      const Vec3d & offs){
    Vec2f co2 = warp(x, y, disp, offs);
    Vec2f co = warp_inv(co2.x, co2.y, flow.bilinear(x, y, 0), offs*-1);
    
    return co - Vec2f(x, y);
  }
  
  /**
     \brief This constraint measures forward-backward flow.
            
     co1 = P_1(bp_0(x, y, d, offs))
     offs1 = Offs1(co1.x, co1.y)
     disp1 = Disp1(co1.x, co1.y)
     
     returns P_0(bp_1(co1.x, co1.y, disp1)) - [x;y]

     // in code notation
     P_0(bp_1(co1.x, co1.y, disp1)) = warp_inv(co1.x, co1.y, disp1, offs1)

     co1 = warp(co.x, co.y, disp, offs)

     partial_co1

  */
  Vec2f forwardBackward(const DisparityOffsetBase<UseDisparity> & warp,
			const DisparityOffsetBase<UseDisparity> & warp_inv,
			const Imagef & flow,
			int x, int y,
			double disp,
			const Vec3d & offs){
    Vec2f co2 = warp(x, y, disp, offs);
    Vec3f offs2(flow.bilinear(co2.x, co2.y, 1), 
		flow.bilinear(co2.x, co2.y, 2),
		flow.bilinear(co2.x, co2.y, 3));
		
    Vec2f co = warp_inv(co2.x, co2.y, flow.bilinear(co2.x, co2.y, 0), offs2);    
    return co - Vec2f(x, y);
  }

  /**
     Quick test: disparity is enforced to not move.
  */
  void addDisparityConstraint(float _alpha = 100.0){
    int S_base = S_I.size();
    int width = f_h.w, height = f_h.h;
    
    S_I.push_back(nacb::Imagef(width, height, 15));
    S_alpha.push_back(_alpha);
    
    S_I.back() = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	Vec5d II(1, 0, 0, 0, 0);
	
	Mat5x5 mat = Mat5x5::outerProduct(II);
		
	if(!mat.isNormal())
	  continue;
		
	accumulate(S_I[S_base], x, y, Mat5x5::outerProduct(II));	
      }
    } 
  }

  Vec2d simpleFlowConstraint2D(int x, int y, 
			       double d, const Vec3d & offs,
			       const nacb::Vec2d & co, 
			       const nacb::Vec2d co_deriv[4], 
			       const DisparityOffsetBase<UseDisparity> & warp, 
			       const DisparityOffsetBase<UseDisparity> & warp_inv,
			       const nacb::Imagef & flow,
			       const Vec4f & Ix_co,
			       const Vec4f & Iy_co,
			       nacb::Vec2d derivs[4]){
    nacb::Vec2d co2_simple_deriv[4];
    nacb::Vec2d co2_simple_xy[2];

    float codisp = flow.bilinear(co.x, co.y, 0);
    Vec2d co2 =  warp_inv(co.x, co.y, codisp, offs*-1, co2_simple_deriv, co2_simple_xy);
    Vec2d result = co2 - nacb::Vec2d(x, y);

    if(derivs && co_deriv){
      for(int j=0; j<4; j++){
	// Derivative of the simple cost function.
	derivs[j] = co2_simple_xy[0]*co_deriv[j].x + co2_simple_xy[1]*co_deriv[j].y +
	  co2_simple_deriv[0]*(Ix_co.data[0]*co_deriv[j].x +
			       Iy_co.data[0]*co_deriv[j].y);
	
	// Non of the inverse_warp coordinates are a function of the variables (hence simpler)
	if(j>=1)
	  derivs[j] -= co2_simple_deriv[j];
      }

#ifdef FB_CHECK_DERIVS
      double eps = 1e-4;
      nacb::Vec2d f_upd[4];
      std::cout << __FUNCTION__ << " deriv check\n";
      for(int j=0; j<4; j++) {
	double dup = d + (j==0?eps:0);
	Vec3d offsup = offs + Vec3d(j==1, j==2, j==3)*eps;
	f_upd[j] = (simpleFlowConstraint2D(x, y, dup, offsup, warp(x, y, dup, offsup), 0, 
					   warp, warp_inv, flow, Ix_co, Iy_co, 0) - result)*(1.0/eps);

	std::cout << derivs[j] << ", " << f_upd[j] << "\n";
      }
#endif
    }

    return result;
  }

  Vec2d forwardBackwardConstraint2D(int x, int y, 
				    double d, const Vec3d & offs,
				    const nacb::Vec2d & co, 
				    const nacb::Vec2d * co_deriv, 
				    const DisparityOffsetBase<UseDisparity> & warp, 
				    const DisparityOffsetBase<UseDisparity> & warp_inv,
				    const nacb::Imagef & flow,
				    const Vec4f & Ix_co,
				    const Vec4f & Iy_co,
				    nacb::Vec2d * derivs){
    float codisp = flow.bilinear(co.x, co.y, 0);
    nacb::Vec2d co2_deriv[4];
    nacb::Vec2d co2_xy[2];
    nacb::Vec3d offs2(flow.bilinear(co.x, co.y, 1),
		      flow.bilinear(co.x, co.y, 2),
		      flow.bilinear(co.x, co.y, 3));
    nacb::Vec2d co2 = warp_inv(co.x, co.y, codisp, offs2, co2_deriv, co2_xy);
    nacb::Vec2d result = co2 - nacb::Vec2d(x, y);

    if(derivs && co_deriv){
      for(int j=0; j<4; j++){
	derivs[j] = co2_xy[0]*co_deriv[j].x + co2_xy[1]*co_deriv[j].y +
	  co2_deriv[0]*(Ix_co.data[0]*co_deriv[j].x +
			Iy_co.data[0]*co_deriv[j].y) +
	  co2_deriv[1]*(Ix_co.data[1]*co_deriv[j].x +
			Iy_co.data[1]*co_deriv[j].y) +
	  co2_deriv[2]*(Ix_co.data[2]*co_deriv[j].x +
			Iy_co.data[2]*co_deriv[j].y) +
	  co2_deriv[3]*(Ix_co.data[3]*co_deriv[j].x +
			Iy_co.data[3]*co_deriv[j].y);      
      }

#ifdef FB_CHECK_DERIVS
      double eps = 1e-4;
      nacb::Vec2d f_upd[4];
      for(int j=0; j<4; j++) {
	double dup = d + (j==0?eps:0);
	Vec3d offsup = offs + Vec3d(j==1, j==2, j==3)*eps;
	f_upd[j] = (forwardBackwardConstraint2D(x, y, dup, offsup, warp(x, y, dup, offsup), 0, 
						warp, warp_inv, flow, Ix_co, Iy_co, 0) - result)*(1.0/eps);

	std::cout << derivs[j] << ", " << f_upd[j] << "\n";
      }
#endif
    }
    return result;
  }


  Vec3d forwardBackwardConstraint3D(int x, int y, 
				    double d, const Vec3d & offs,
				    const nacb::Vec2d & co, 
				    const nacb::Vec2d * co_deriv, 
				    const DisparityOffsetBase<UseDisparity> & warp, 
				    const DisparityOffsetBase<UseDisparity> & warp_inv,
				    const nacb::Imagef & flow,
				    const Vec4f & Ix_co,
				    const Vec4f & Iy_co,
				    nacb::Vec3d * derivs){
    float codisp = flow.bilinear(co.x, co.y, 0);
    nacb::Vec3d co2_deriv[4];
    nacb::Vec3d co2_xy[2];
    nacb::Vec3d offs2(flow.bilinear(co.x, co.y, 1),
		      flow.bilinear(co.x, co.y, 2),
		      flow.bilinear(co.x, co.y, 3));
    nacb::Vec3d backProjDeriv;
    nacb::Vec3d co2 = warp_inv.backProjectAndOffset(co.x, co.y, codisp, offs2, co2_deriv, co2_xy);
    nacb::Vec3d result = co2 - warp.backProject(x, y, d, &backProjDeriv);

    if(derivs && co_deriv){
      for(int j=0; j<4; j++){
	derivs[j] = co2_xy[0]*co_deriv[j].x + co2_xy[1]*co_deriv[j].y +
	  co2_deriv[0]*(Ix_co.data[0]*co_deriv[j].x +
			Iy_co.data[0]*co_deriv[j].y) +
	  co2_deriv[1]*(Ix_co.data[1]*co_deriv[j].x +
			Iy_co.data[1]*co_deriv[j].y) +
	  co2_deriv[2]*(Ix_co.data[2]*co_deriv[j].x +
			Iy_co.data[2]*co_deriv[j].y) +
	  co2_deriv[3]*(Ix_co.data[3]*co_deriv[j].x +
			Iy_co.data[3]*co_deriv[j].y);      
      }
      derivs[0] -= backProjDeriv;

#ifdef FB_CHECK_DERIVS_3D
      double eps = 1e-4;
      nacb::Vec3d f_upd[4];

      std::cout << "fb3D deriv check\n";

      for(int j=0; j<4; j++) {
	double dup = d + (j==0?eps:0);
	Vec3d offsup = offs + Vec3d(j==1, j==2, j==3)*eps;
	
	double ddn = d - (j==0?eps:0);
	Vec3d offsdn = offs - Vec3d(j==1, j==2, j==3)*eps;

	f_upd[j] = (forwardBackwardConstraint3D(x, y, dup, offsup, warp(x, y, dup, offsup), 0, 
						warp, warp_inv, flow, Ix_co, Iy_co, 0) - 
		    forwardBackwardConstraint3D(x, y, ddn, offsdn, warp(x, y, ddn, offsdn), 0, 
						warp, warp_inv, flow, Ix_co, Iy_co, 0))*(0.5/eps);

	std::cout << derivs[j] << ", " << f_upd[j] << "\n";
      }
#endif
    }
    return result;
  }


  /**
     Works okay.  Not great.
   */
  void addDisplacementConsistency(const std::vector<DisparityOffsetBase<UseDisparity> > & warps,
				  const std::vector<DisparityOffsetBase<UseDisparity> > & warp_invs,
				  const std::vector<Imagef> & flows,
				  float _alpha = 1.0){
    bool constrain2D = false;
    int S_base = S_I.size();

    // Add the extra terms.
    int width = f_h.w, height = f_h.h;

    for(int i=0; i<(int)warps.size(); i++){
      S_I.push_back(nacb::Imagef(width, height, 15));
      S_alpha.push_back(_alpha);
      
      S_I.back() = 0;
    }

    nacb::Imagef res(width, height, 1);
    res = 0;

    double totalDiff = 0;
    int ndiff = 0;

    // Compute the actual derivatives.
    for(int i=0; i<(int)warps.size(); i++){
      const DisparityOffsetBase<UseDisparity> & warp = warps[i];
      const nacb::Imagef & flow = flows[i];

      nacb::Imagef Ix, Iy;
      image_gradient_all(flow, Ix, Iy);
      Ix *= (1.0/hx);
      Iy *= (1.0/hy);
      
      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  if(weight(x, y) >= 1.0){
	    // Compute the offset coordinates, and its derivatives.
	    nacb::Vec2d co, co_deriv[4];
	    nacb::Vec3d offs(D_0(x, y, 1), D_0(x, y, 2), D_0(x, y, 3));
	    co = warp(x, y, D_0(x, y, 0), offs, co_deriv);

	    // Sample the neighbors flow at the offset coords
	    nacb::Vec4f Ix_co, Iy_co;
	    Ix.bilinear(co.x, co.y, Ix_co.data, 4);
	    Iy.bilinear(co.x, co.y, Iy_co.data, 4);

	    if(constrain2D){
	      // Add the 2D consistency constraints.
	      for(int funci=0; funci<2; funci++){
		Vec2d codiff;
		Vec2d dvs[4];
		
		if(funci == 0)
		  codiff = forwardBackwardConstraint2D(x, y, D_0(x, y, 0), offs, co, co_deriv, 
						       warps[i], warp_invs[i], flow, Ix_co, Iy_co, dvs);
		else 
		  codiff = simpleFlowConstraint2D(x, y, D_0(x, y, 0), offs, co, co_deriv, 
						  warps[i], warp_invs[i], flow, Ix_co, Iy_co, dvs);
		
		for(int k=0; k<2; k++){
		  double Iz = codiff.data[k];
		  Vec5d II(dvs[0].data[k], dvs[1].data[k],
			   dvs[2].data[k], dvs[3].data[k], Iz);
		  
		  Mat5x5 mat = Mat5x5::outerProduct(II);
		  
		  if(!mat.isNormal())
		    continue;
		  
		  accumulate(S_I[i + S_base], x, y, Mat5x5::outerProduct(II));
		  
		  res(x, y) += Iz*Iz;
		  totalDiff += weight(x, y)*Iz*Iz;
		  ndiff++;
		}
	      }
	    }
	    else {
	      // 3D consisntency

	      Vec3d dvs[4];
	      Vec3d codiff = forwardBackwardConstraint3D(x, y, D_0(x, y, 0), offs, co, co_deriv, 
							 warps[i], warp_invs[i], flow, Ix_co, Iy_co, dvs);
	      
	      for(int k=0; k<3; k++){
		double Iz = codiff.data[k];
		Vec5d II(dvs[0].data[k], dvs[1].data[k],
			 dvs[2].data[k], dvs[3].data[k], Iz);
	      
		Mat5x5 mat = Mat5x5::outerProduct(II);
		
		if(!mat.isNormal())
		  continue;
		
		accumulate(S_I[i + S_base], x, y, Mat5x5::outerProduct(II));
		
		res(x, y) += Iz*Iz;
		totalDiff += weight(x, y)*Iz*Iz;
		ndiff++;
	      }
	      
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


  VariationalSceneFlowProblem  restriction();
  
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
  static double imageDataEnergy(const Imagef & D,
                                const nacb::Imagef & image,
                                const nacb::Imagef & weight,
				const std::vector<DisparityOffsetBase<UseDisparity> > & warps, 
                                const std::vector<nacb::Imagef> & im2s,
                                double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
        
	nacb::Vec3f color(image(x, y, 0), image(x, y, 1), image(x, y, 2));
        
        for(int si=0; si<(int)warps.size(); si++){
	  Vec3f offs(D(x, y, 1), D(x, y, 2), D(x, y, 3));
          Vec2f co = warps[si](x, y, D(x, y, 0), offs);
          Vec3f color2;
          
	  im2s[si].bilinear(co.x, co.y, color2.data, 3);
          
          Vec3f diff = color - color2;
          energy_D += weight(x, y)*phi_D(diff.dot(diff));
        }       
      }
    }
    return energy_D;
  }


  /**
     Compute the geometric data energy (3D)
  */
  static double geometricDataEnergy2D(const Imagef & D,
				      const nacb::Imagef & image,
				      const nacb::Imagef & weight,
				      const std::vector<DisparityOffsetBase<UseDisparity> > & warps, 
				      const std::vector<DisparityOffsetBase<UseDisparity> > & warp_invs,
				      const std::vector<nacb::Imagef> & im2s,
				      const std::vector<nacb::Imagef> & flows,
				      double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
        
	nacb::Vec3f color(image(x, y, 0), image(x, y, 1), image(x, y, 2));
        
        for(int si=0; si<(int)warps.size(); si++){
	  Vec3f offs(D(x, y, 1), D(x, y, 2), D(x, y, 3));
          Vec2f co = warps[si](x, y, D(x, y, 0), offs);
	  Vec4f flowOffs;

	  flows[si].bilinear(co.x, co.y, flowOffs.data, 4);
	  
	  Vec2f cons_simple = 
	    warp_invs[si](co.x, co.y, flowOffs.x, offs*-1) - Vec2f(x, y);
	  
	  Vec2f cons_sym = 
	    warp_invs[si](co.x, co.y, flowOffs.x, Vec3f(flowOffs.y, flowOffs.z, flowOffs.w)) - Vec2f(x, y);

	  energy_D += weight(x, y)*phi_D(cons_simple.dot(cons_simple));
	  energy_D += weight(x, y)*phi_D(cons_sym.dot(cons_sym));
        }       
      }
    }
    return energy_D;
  }


  /**
     Compute the geometric data energy (3D)
  */
  static double geometricDataEnergy3D(const Imagef & D,
				      const nacb::Imagef & image,
				      const nacb::Imagef & weight,
				      const std::vector<DisparityOffsetBase<UseDisparity> > & warps, 
				      const std::vector<DisparityOffsetBase<UseDisparity> > & warp_invs,
				      const std::vector<nacb::Imagef> & im2s,
				      const std::vector<nacb::Imagef> & flows,
				      double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
        
	nacb::Vec3f color(image(x, y, 0), image(x, y, 1), image(x, y, 2));
        
        for(int si=0; si<(int)warps.size(); si++){
	  Vec3f offs(D(x, y, 1), D(x, y, 2), D(x, y, 3));
          Vec2f co = warps[si](x, y, D(x, y, 0), offs);
	  Vec4f flowOffs;

	  
	  flows[si].bilinear(co.x, co.y, flowOffs.data, 4);
	  
	  Vec3d bpoffs2 = warp_invs[si].backProjectAndOffset(co.x, co.y, flowOffs.x, 
							     Vec3f(flowOffs.y, flowOffs.z, flowOffs.w));
	  Vec3d bp = warps[si].backProject(x, y, D(x, y, 0));
	  Vec3d cons = bpoffs2 - bp;

	  energy_D += weight(x, y)*phi_D(cons.dot(cons));
        }       
      }
    }
    return energy_D;
  }
};



#endif // VARIATIONAL_SF_H
