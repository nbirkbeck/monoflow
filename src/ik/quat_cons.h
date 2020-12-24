#ifndef QUAT_CONS_H
#define QUAT_CONS_H

#include <assert.h>
#include <nmath/matrix.h>
#include <nmath/quaternion.h>
#include <math.h>

double boundedClosestXYZ(const nacb::Quaternion & q,
			 double x[],
			 double l[],
			 double u[]);
#endif
