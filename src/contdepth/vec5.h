#ifndef NMATH_VEC5_H
#define NMATH_VEC5_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

namespace nacb {

  class Vec5d {

  public:
    Vec5d(double i = 0, double j = 0, double k = 0, double l = 0, double m = 0) {
      data[0] = i;
      data[1] = j;
      data[2] = k;
      data[3] = l;
      data[4] = m;
    }

    Vec5d operator*(double w) const {
      Vec5d res;
      for(int j=0; j<5; j++)
	res.data[j] = data[j]*w;
      return res;
    }

    Vec5d & operator+=(const Vec5d & o) {
      for(int i=0; i<5; i++)
	data[i] += o.data[i];
      return *this;
    }

    double & operator[](int j) {
      return data[j];
    }

    const double & operator[](int j) const {
      return data[j];
    }

    double dot(const Vec5d & o) const {
      double sum = 0;
      for(int i=0; i<5; i++)
	sum += o.data[i]*data[i];
      return sum;
    }

    union {
      double data[5];
      struct {
	double x, y, z, w, a;
      };
    };
  };
  

  class Mat5x5 {

  public:
    Mat5x5() { }

    Vec5d operator*(const Vec5d & other) const {
      Vec5d res;
      
      for(int i=0; i<5; i++){
	for(int j=0; j<5; j++)
	  res.data[i] += data[i][j]*other.data[j];
      }

      return res;
    }

    bool isNormal() const {
      for(int ii=0; ii<5; ii++){
	for(int jj=0; jj<5; jj++){
	  if(isnan(data[ii][jj]) || isinf(data[ii][jj]))
	    return false;
	}
      }
      return true;
    }

    double & operator()(int i, int j) {
      return data[i][j];
    }

    const double & operator()(int i, int j) const {
      return data[i][j];
    }

    static Mat5x5 zeros() {
      Mat5x5 m;
      memset(m.data, 0, sizeof(m));
      return m;
    }

    static Mat5x5 outerProduct(const Vec5d & v) {
      Mat5x5 mat = Mat5x5::zeros();
     
      for(int i=0; i<5; i++)
	for(int j=0; j<5; j++)
	  mat(i, j) += v[i]*v[j];

      return mat;
    }

    double data[5][5];
  };
};

#endif //NMATH_VEC5_H
