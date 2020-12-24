#ifndef HERMITE_H
#define HERMITE_H


double hermite_tangent(const float * data, int i0, int n);

double hermite_interp(const float * data, float index, 
		      int n, bool deriv = false);

#endif // HERMITE_H
