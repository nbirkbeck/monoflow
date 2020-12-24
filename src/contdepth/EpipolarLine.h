#ifndef EPIPOLAR_LINE_H
#define EPIPOLAR_LINE_H

#include <nmath/matrix.h>
#include <nmath/vec2.h>
#include <nmath/vec3.h>

template <bool useDisparity = true>
struct EpipolarLineTemplate {
  nacb::Matrix A;
  nacb::Matrix b;

  EpipolarLineTemplate(const nacb::Matrix & K1, const nacb::Matrix & E1,
		       const nacb::Matrix & K2, const nacb::Matrix & E2){
    A = K2*E2.submatrix(0, 0, 3, 3)*E1.submatrix(0, 0, 3, 3).transpose()*K1.inverse();

    nacb::Matrix E1inv = E1.inverse();
    nacb::Matrix E2inv = E2.inverse();

    b = K2*E2.submatrix(0, 0, 3, 3)*(E1inv.submatrix(0, 3, 3, 1) - E2inv.submatrix(0, 3, 3, 1));
  }

  nacb::Vec3<float> mult(float x, float y, float d) const {
    if (useDisparity) {
      return nacb::Vec3<float>(A(0, 0)*x + A(0, 1)*y + A(0, 2) + d*b[0],
			       A(1, 0)*x + A(1, 1)*y + A(1, 2) + d*b[1],
			       A(2, 0)*x + A(2, 1)*y + A(2, 2) + d*b[2]);
    }
    else {
      return nacb::Vec3<float>(A(0, 0)*x + A(0, 1)*y + A(0, 2) + b[0]/d,
			       A(1, 0)*x + A(1, 1)*y + A(1, 2) + b[1]/d,
			       A(2, 0)*x + A(2, 1)*y + A(2, 2) + b[2]/d);
    }

  }

  nacb::Vec2<float> operator()(float x, float y, float d) const {
    nacb::Vec3<float> p = mult(x, y, d);
    return nacb::Vec2<float>(p.x/p.z, p.y/p.z);
  }

  nacb::Vec2<float> disparityDeriv(float x, float y, float d) const {
    nacb::Vec3<float> co = mult(x, y, d);
    if (useDisparity) {
      return nacb::Vec2<float>(b[0]/co.z - co.x*b[2]/(co.z*co.z),
			       b[1]/co.z - co.y*b[2]/(co.z*co.z));
    }
    else {
      return nacb::Vec2<float>(-b[0]/(co.z*d*d) + co.x*b[2]/(co.z*co.z*d*d),
			       -b[1]/(co.z*d*d) + co.y*b[2]/(co.z*co.z*d*d));
    }
  }

  void deriv_xy(float x, float y, float d, nacb::Vec2f * dxy) const {
    nacb::Vec3<float> co = mult(x, y, d);
    dxy[0] = nacb::Vec2f(A(0, 0)/co.z - A(2, 0)*co.x/(co.z*co.z),
			 A(0, 1)/co.z - A(2, 1)*co.x/(co.z*co.z));
    
    dxy[1] = nacb::Vec2f(A(1, 0)/co.z - A(2, 0)*co.y/(co.z*co.z),
			 A(1, 1)/co.z - A(2, 1)*co.y/(co.z*co.z));
    
  }

  static const bool UseDisparity = useDisparity; 
};


typedef EpipolarLineTemplate<true> EpipolarLine;

#endif //EPIPOLAR_LINE_H
