#ifndef DISPARITY_OFFSET_H
#define DISPARITY_OFFSET_H

#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"

#include <nmath/matrix.h>
#include <nmath/vec2.h>
#include <nmath/vec3.h>

/**
   \brief DisparityOffsetBase is a class to hold the displacement
          warp (for a disparity and 3D displacement).

	  useDisparity: true if we are using disparity
	                otherwise we are using depth-based formulation.
*/
template <bool useDisparity = true>
class DisparityOffsetBase {
 public:

  DisparityOffsetBase(const nacb::Matrix & A1, const nacb::Matrix & E1,
		      const nacb::Matrix & A2, const nacb::Matrix & E2,
		      double _offsScale = 1): offsetScale(_offsScale) {
    nacb::Matrix A = nacb::Matrix::eye(4, 4);
    A.setSubmatrix(0, 0, A1.inverse());
    E1i = E1.inverse() * A;
    P = A2 * nacb::Matrix::eye(3, 4) * E2;
    PE1i = P*E1i;
  }

  nacb::Vec3d backProject(float x, float y, float disp, 
			  nacb::Vec3d * derivs = 0,
			  nacb::Vec3d * dxys = 0) const {
    nacb::Matrix X(4, 1);

    if (useDisparity) {
      X[0] = x/disp;
      X[1] = y/disp;
      X[2] = 1.0/disp;
      X[3] = 1;
    }
    else {
      X[0] = x*disp;
      X[1] = y*disp;
      X[2] = disp;
      X[3] = 1;
    }

    X = E1i * X;

    if(derivs){
      nacb::Matrix ddisp(4, 1);

      if (useDisparity) {
	ddisp[0] = -x/(disp*disp);
	ddisp[1] = -y/(disp*disp);
	ddisp[2] = -1.0/(disp*disp);
	ddisp[3] = 0;
      }
      else {
	ddisp[0] = x;
	ddisp[1] = y;
	ddisp[2] = 1;
	ddisp[3] = 0;
      }

      nacb::Matrix E1i_ddisp = E1i*ddisp;

      derivs[0] = nacb::Vec3d(E1i_ddisp[0], E1i_ddisp[1], E1i_ddisp[2]);
    }
    if(dxys){
      if (useDisparity) {
	dxys[0] = nacb::Vec3d(E1i(0, 0)/disp, E1i(1, 0)/disp, E1i(2, 0)/disp);
	dxys[1] = nacb::Vec3d(E1i(0, 1)/disp, E1i(1, 1)/disp, E1i(2, 1)/disp);
      }
      else {
	dxys[0] = nacb::Vec3d(E1i(0, 0)*disp, E1i(1, 0)*disp, E1i(2, 0)*disp);
	dxys[1] = nacb::Vec3d(E1i(0, 1)*disp, E1i(1, 1)*disp, E1i(2, 1)*disp);
      }
    }
    
    return nacb::Vec3d(X[0], X[1], X[2]);
  }

  nacb::Vec3d backProjectAndOffset(float x, float y, float disp, const nacb::Vec3d & offs,
				   nacb::Vec3d * derivs = 0,
				   nacb::Vec3d * dxys = 0) const {
    nacb::Vec3d bp = backProject(x, y, disp, derivs, dxys);
    
    // Back-project is not a function of offsets, and the
    // rest of this function is not a function of x, y or disp
    // Only update the offset derivaties.
    if(derivs){
      derivs[1] = nacb::Vec3d(offsetScale, 0, 0);
      derivs[2] = nacb::Vec3d(0, offsetScale, 0);
      derivs[3] = nacb::Vec3d(0, 0, offsetScale);
    }
    return offs*offsetScale + bp;
  }


  nacb::Vec2d operator()(float x, float y, float disp, 
		   const nacb::Vec3d & offs, 
		   nacb::Vec2d * derivs = 0,
		   nacb::Vec2d * dxys = 0) const {
    /*
    nacb::Matrix X(4, 1);

    X[0] = x/disp;
    X[1] = y/disp;
    X[2] = 1.0/disp;
    X[3] = 1;

    X = E1i * X;
    
    nacb::Vec3d co = offs*offsetScale + nacb::Vec3d(X[0], X[1], X[2]);
    */
    nacb::Vec3d co = backProjectAndOffset(x, y, disp, offs);
    double pz = 0;
    nacb::Vec2d pj = project(P, co, &pz);

    // Compute the derivative if the derivs is here.
    if(derivs){
      nacb::Vec3d xhat(pj.x*pz, pj.y*pz, pz);
      double z2 = pz*pz;

      nacb::Vec3d dxhat[4];

      nacb::Matrix ddisp(4, 1);

      if (useDisparity) {
	ddisp[0] = -x/(disp*disp);
	ddisp[1] = -y/(disp*disp);
	ddisp[2] = -1.0/(disp*disp);
	ddisp[3] = 0;
      }
      else {
	ddisp[0] = x;
	ddisp[1] = y;
	ddisp[2] = 1.0;
	ddisp[3] = 0;
      }

      nacb::Matrix PE1i_ddisp = PE1i*ddisp;
      
      dxhat[0] = nacb::Vec3d(PE1i_ddisp[0], PE1i_ddisp[1], PE1i_ddisp[2]);
      dxhat[1] = nacb::Vec3d(P(0, 0), P(1, 0), P(2, 0))*offsetScale;
      dxhat[2] = nacb::Vec3d(P(0, 1), P(1, 1), P(2, 1))*offsetScale;
      dxhat[3] = nacb::Vec3d(P(0, 2), P(1, 2), P(2, 2))*offsetScale;
      
      for(int j=0; j<4; j++)
	derivs[j] = nacb::Vec2d((xhat.z * dxhat[j].x - xhat.x * dxhat[j].z)/z2,
			  (xhat.z * dxhat[j].y - xhat.y * dxhat[j].z)/z2);

      // Partial derivative of the x coord w.r.t the x,y
      if(dxys){
	for(int j=0; j<2; j++){
	  if (useDisparity) {
	    nacb::Vec3d dxhat = nacb::Vec3d(PE1i(0, j), PE1i(1, j), PE1i(2, j))*(1.0/disp);
	    
	    dxys[j] = nacb::Vec2d((xhat.z * dxhat.x - xhat.x * dxhat.z)/z2,
				  (xhat.z * dxhat.y - xhat.y * dxhat.z)/z2);
	  }
	  else {
	    nacb::Vec3d dxhat = nacb::Vec3d(PE1i(0, j), PE1i(1, j), PE1i(2, j))*disp;
	    
	    dxys[j] = nacb::Vec2d((xhat.z * dxhat.x - xhat.x * dxhat.z)/z2,
				  (xhat.z * dxhat.y - xhat.y * dxhat.z)/z2);
	  }
	}
      }

#ifdef CHECK_DERIV
      double eps = 1e-4;
      nacb::Vec2d ups[4] = {operator()(x, y, disp + eps, offs),
			    operator()(x, y, disp, offs + nacb::Vec3d(eps, 0, 0)),
			    operator()(x, y, disp, offs + nacb::Vec3d(0, eps, 0)),
			    operator()(x, y, disp, offs + nacb::Vec3d(0, 0, eps))};
      for(int j=0; j<4; j++){
	ups[j] -= pj;
	ups[j] *= (1.0/eps);
	
	printf("%d %f %f  %f %f  %f\n", j,
	       derivs[j].x, ups[j].x, derivs[j].y, ups[j].y, 
	       (ups[j] - derivs[j]).len());
      }
#endif
    }

    return pj;
  }

  static const bool UseDisparity = useDisparity; 
  double offsetScale;
  nacb::Matrix E1i, P, PE1i;
};

typedef DisparityOffsetBase<true> DisparityOffset;

#endif // DISPARITY_OFFSET
