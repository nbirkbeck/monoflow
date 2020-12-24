#ifndef RECON_GEOMETRY_H
#define RECON_GEOMETRY_H

#include "feature.h"
#include <nimage/image.h>
#include <nmath/matrix.h>
#include <vector>

using nacb::Matrix;
using nacb::Image8;

Matrix skewSymmetric(const Matrix & v);

Matrix triangulate(vector<Matrix> & Ps,
		   vector<Matrix> & xs);

double triangulate_opt(vector<Matrix> & Ps,
		       vector<Matrix> & xs,
		       Matrix & X);

template <class T>
nacb::Vec3<T> triangulate(const nacb::Matrix & P1, const nacb::Matrix & P2,
			  const nacb::Vec2<T> & p1, const nacb::Vec2<T> & p2);

std::vector<int>  triangulate_robust(std::vector<nacb::Matrix> & Ps,
				     std::vector<nacb::Matrix> & xs,
				     nacb::Matrix & X);

vector<int> alignRobust(Matrix & H,vector<Matrix> & Ps,
			Matrix & X,Matrix & xs,vector<int> & cams,
			int euclidean=0);

void align_refine(Matrix & H,vector<Matrix> & Ps,
		  Matrix & X,Matrix & xs,vector<int> & cams);

Matrix getFundamentalMatrix(Matrix & p1in,Matrix & p2in);
vector<int> fmatrix(std::vector<Feature *> & f1,	     
		    std::vector<Feature *> & f2,
		    double ethresh = 1.5,
		    double distanceLimit= 100000);

//returns projection matrices
vector<Matrix> robustTrilinear(Matrix & ps1,Matrix & ps2,Matrix & ps3,
			       Matrix * Xout=0,
			       vector<int> * inds_out=0);
vector<Matrix>  getTrilinearTensorMinimal(Matrix & p1in,Matrix & p2in,Matrix & p3in,
					  Matrix & p1eval,//all the points to evaluate the 3 soln's
					  Matrix & p2eval,
					  Matrix & p3eval,
					  int * bestscore_out=0);

double homographyModelError(Matrix & x1, Matrix & x2);

Matrix getHomography(Matrix & X_,Matrix & XP_,int dlt=0);

Matrix projectiveFactorization(vector<Matrix> & obs,
			       vector<Matrix> & Ps);

Matrix resect(Matrix & X, Matrix & x);
Matrix resectEuclidean(Matrix & K,Matrix & X,Matrix & x, bool dofocal=false);


//From svn/thesis/calibration
int calibratePlanar(Matrix pattern,
		    Matrix * pts,int nimages,
		    Matrix & A,
		    Matrix & dist,
		    Matrix * Rs,
		    Matrix * ts);

#endif
