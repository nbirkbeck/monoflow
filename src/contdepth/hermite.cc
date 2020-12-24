#include "hermite.h"
#include <nmath/matrix.h>
#include <stdio.h>

using namespace nacb;

double hermite_tangent(const float * data, int i0, int n){
  if(i0 == 0)
    return data[1] - data[0];
  if(i0 == (n - 1))
    return data[n - 1] - data[n - 2];

  return (data[i0 + 1] - data[i0 - 1])/2;
}


double hermite_interp(const float * data, float index, 
		     int n, bool deriv){
  if(index < 0){
    if(deriv)
      return -1;
    return data[0] - index;
  }
  else if (index >= n){
    if(deriv)
      return 1;
    return (index - (n - 1)) + data[n - 1];
  }
  
  int i0 = (int)floor(index);
  int i1 = i0 + 1;

  double t = index - i0;
  double m0 = hermite_tangent(data, i0, n);
  double m1 = hermite_tangent(data, i1, n);

  double a, b, c, d;
  if(deriv){
    a = 6.*t*t - 6.*t;
    b = 3.*t*t - 4.*t + 1.;
    c = -6.*t*t + 6.*t;
    d = 3.*t*t - 2.*t;
  }
  else {
    a = 2.*t*t*t - 3.*t*t + 1.;
    b = t*t*t - 2.*t*t + t;
    c = -2.*t*t*t + 3.*t*t;
    d = t*t*t - t*t;
  }

  return data[i0]*a + m0*b + data[i1]*c + m1*d;
}



void test_hermite_interp(){
  float x[] = {1.000000,1.050000,1.100000, 1.150000,1.200000,1.250000,
	       1.300000,1.350000,1.400000, 1.450000,1.500000,1.550000,
	       1.600000,1.650000,1.700000, 1.750000,1.800000,1.850000,
	       1.900000,1.950000,2.000000};
  float c[] = {2.1, 2, 1.9, 2.2}; 
  Matrix y(21, 1);

  for(int i=0; i<21; i++){
    y[i] = hermite_interp(c, x[i], 4);
    printf("grad: %f %f\n", hermite_interp(c, x[i], 4, true), 
	   (hermite_interp(c, x[i] + 1e-5, 4) - y[i])/1e-5);
  }
  y.printMatlab("y");
}
