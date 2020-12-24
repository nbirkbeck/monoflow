#ifndef BUNDLE_ADJUST_H
#define BUNDLE_ADJUST_H

#include <nmath/matrix.h>
#include "feature.h"
#include <vector>
#include <string>

using nacb::Matrix;

double bundleAdjustSparse(std::vector<Matrix> & Ps,
			  std::vector<Matrix> & xs,
			  Matrix & X);

void bundleAdjustSparse(std::vector<Matrix> & Ps,
			std::vector<std::vector<Feature *> > & features,
			std::vector<std::vector<int> > & x2X,
			std::vector<Matrix> & X);
//euclidean BA
void bundleAdjustSparse(std::vector<Matrix> & Ks,
			std::vector<Matrix> & Rts,
			std::vector<std::vector<Feature *> > & features,
			std::vector<std::vector<int> > & x2X,
			Matrix & X,int dofocal=1);

//no need to be sparse for one camera
double bundleAdjustCamera(const Matrix & K,
			  Matrix & Rt,
			  const Matrix & xs,
			  const Matrix & X,int dofocal);
//no need to be sparse for one camera
double bundleAdjustCamera(const Matrix & K,
			Matrix & Rt,
			std::vector<Feature *> & features,
			std::vector<int> & x2X,
			const Matrix & X,int dofocal=1);
void bundleAdjust(std::vector<Matrix> & Ps,
		  std::vector<std::vector<Feature *> > & features,
		  std::vector<std::vector<int> > & x2X,
		  std::vector<Matrix> & X);
double bundleAdjust(std::vector<Matrix> & Ps,
		    std::vector<Matrix> & xs,
		    Matrix & X);
#endif
