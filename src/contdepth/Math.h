#ifndef CONTDEPTH_MATH_H
#define CONTDEPTH_MATH_H

#include <math.h>
#include "vec5.h"
#include <nimage/image.h>
#include <nmath/matrix.h>
#include <assert.h>
#include <stdexcept>

using namespace nacb;

inline double sqr(double d){
  return d*d;
}

bool isNormal(nacb::Matrix & m);

template <class T>
double dot4(Vec4<T>& t1, Vec4<T>& t2) {
  double sum = 0;
  for (int k=0; k<4; k++)
    sum += t1[k] * t2[k];
  return sum;
}

/**
   \brief Pack a symmetric matrix into a 1D vector.
 */
inline nacb::Matrix packSymmetricMatrix(const nacb::Matrix& minput) {
  nacb::Matrix m((minput.m + 1)*minput.n/2, 1);
  int k = 0;
  for (int i=0; i<minput.m; i++) {
    for (int j=i; j<minput.n; j++, k++) {
      m[k] = minput(i, j);
    }
  }
  return m;
}


/** \brief Unpack a symmetric matrix from a 1D vector into a matrix.
   The reason this doesn't return a matrix is to save on malloc/free.
 */
template <class T>
inline void unpackSymmetricMatrix(const T * mpacked, int npacked, nacb::Matrix& mout) {
  int sz = int(-1 + sqrt(1 + 8*npacked))/2;
  if (mout.m != sz || mout.n != sz)
    mout = nacb::Matrix(sz, sz);

  int k = 0;
  for (int i=0; i<mout.m; i++) {
    for (int j=i; j<mout.n; j++, k++) {
      mout(i, j) = mpacked[k];
      mout(j, i) = mpacked[k];
    }
  }
}


nacb::Mat5x5 unpack_sym_mat5x5(const nacb::Imagef & image, int x, int y);


void solve_symmetric_4x4(Matrix & M, Matrix & b, Matrix & x);
void solve_symmetric_nxn(Matrix & M, Matrix & b, Matrix & x);

inline void matrix_vec_mult(const Matrix& M, const Matrix& x, Matrix& temp) {
 if (temp.m != M.m || temp.n != 1)
    temp = Matrix(M.m, 1);

  // Compute temp = M*x
  for (int ii=0; ii<M.m; ii++) {
    temp[ii] = 0;

    for (int jj=0; jj<M.n; jj++) {
      temp[ii] += M(ii, jj)*x[jj];
    }
  }
}


inline double quadratic_form_eval(const Matrix & M, Matrix & b, Matrix& temp) {
  if (M.m != M.n)
    throw std::runtime_error("Matrix must be symmetric.");

  if (b.n != 1)
    throw std::runtime_error("b must be a column vector");

  if (M.n != b.m)
    throw std::runtime_error("M.m != b.n");

  matrix_vec_mult(M, b, temp);

  // Compute result = b'*M
  return b.dot(temp);
}


inline void copy_matrix_no_alloc(nacb::Matrix & out, const nacb::Matrix & source) {
  for(int i=0; i<out.m*out.n; i++)
    out.data[i] = source.data[i];
}


inline void diff_matrix_no_alloc(nacb::Matrix & out, const nacb::Matrix & a, const nacb::Matrix & b) {
  for(int i=0; i<out.m*out.n; i++)
    out.data[i] = a.data[i] - b.data[i];
}


Imagef get_D_phi_S(const nacb::Imagef U_grad[4], 
		   double hx, double hy,
                   double (*d_phi_S)(double));


/**
   \brief accumulate the matrix into the channels of the symmetric image.
*/
inline void accumulate(nacb::Imagef & S_I, int x, int y, const Mat5x5 & mat){
  int ind = 0;
  
  assert(S_I.nchannels == 15);

  for(int ii=0; ii<5; ii++){
    for(int jj=ii; jj<5; jj++){
      S_I(x, y, ind) += mat(ii, jj);
      ind++;
    }
  }
}


inline void image_gradient_all(const nacb::Imagef & input, nacb::Imagef & gx, nacb::Imagef & gy){
  //#define USE_5_POINT_GRADIENT
#ifdef USE_5_POINT_GRADIENT
  nacb::Imagef fx(5, 1, 1);
  nacb::Imagef fy(1, 5, 1);

  fy(0, 0) = fx(0, 0) =  1.f/12.f;
  fy(0, 1) = fx(1, 0) = -8.f/12.f;
  fy(0, 2) = fx(2, 0) =  0.f;
  fy(0, 3) = fx(3, 0) =  8.f/12.f;
  fy(0, 4) = fx(4, 0) =  1.f/12.f;

  gx = input.convolve(fx);
  gy = input.convolve(fy);

#else
  gx = nacb::Imagef(input.w, input.h, input.nchannels);
  gy = nacb::Imagef(input.w, input.h, input.nchannels);


  for(int k=0; k<input.nchannels; k++){
    nacb::Imagef grad = input.gradient(k);

    gx.setChannel(k, grad.getChannel(0));
    gy.setChannel(k, grad.getChannel(1));
  }
#endif
}

#endif 
