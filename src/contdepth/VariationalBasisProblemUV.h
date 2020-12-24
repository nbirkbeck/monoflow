/**
   // TODO: 
   1) put displacements in the frame of the bone?
    -not sure if this would improve results or not.  I think this is necessary.
    -but, should the displacements be in tangent space??  
      -this probably makes the most sense.
      -the other alternative is to use the bone space, e.g.,
        the linear combination of the transformations.
   2) the basis for the depth does not need to have so many elements.
   3) make regularization epsilon a parameter.

 */
#ifndef VARIATIONAL_H
#define VARIATIONAL_H

#include "DisplaceUV.h"
#include "vec5.h"
#include "Math.h"

#include "autorecon/stereo/flow_variational.h"
#include "autorecon/stereo/multigrid.h"

#include <nmisc/timer.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>

#include <nimage/image.h>
#include <vector>

#include "ik/armature.h"
#include "DisplaceOffsetBasis.h"
#include "FlowBasis3D.h"
#include "BaseMeshAnimation.h"

using namespace nacb;



/** \brief A struct to pack/unpack the non-zero entries of a symmetric matrix
           into a linear array.
    
    For historical reasons it uses the symmetric packing instead of just using the indices
    for both.  This will allow for easier implementations of sparse-matrix multiplication
    (instead of unpacking), if it is necessary.
    
    appendOne is for the homogenous component of the tensor.
 */
struct MatrixPacker {
  int m;
  std::vector<std::pair<int, int> > indices;
  bool dense;
  int numUnknowns;

  MatrixPacker(int numBasis, bool appendOne = true): m(numBasis), dense(true) { 
    if (appendOne) 
      m++;

    numUnknowns = m*(m + 1)/2;
  }
  
  MatrixPacker(const nacb::Matrix& interactVector, bool appendOne = true):dense(false) {
    nacb::Matrix outter;
    if (appendOne) {
      nacb::Matrix M = nacb::Matrix::cat(0, interactVector, nacb::Matrix::ones(1, 1));
      outter = M*M.transpose();
    }
    else
      outter = interactVector*interactVector.transpose();
    
    this->m = outter.m;
    for (int i=0; i<outter.m; i++) {
      for (int j=i; j<outter.n; j++) {
	if (fabs(outter(i, j)) > 1e-8) {
	  indices.push_back(make_pair(i, j));
	}
      }
    }
    numUnknowns = indices.size();
  }

  int getPackedSize() const {
    return numUnknowns;
  }

  nacb::Matrix pack(const nacb::Matrix& M) {
    if (dense) 
      return packSymmetricMatrix(M);
    else {
      nacb::Matrix X(indices.size(), 1);
      for (int i=0; i<(int)indices.size(); i++)
	X[i] = M(indices[i].first, indices[i].second);
      return X;
    }
  }

  void unpack(const float* data, int dataLength, nacb::Matrix& M){
    if (dense) 
      unpackSymmetricMatrix<float>(data, dataLength, M);
    else {
      if (M.m != m || M.n != m) M = nacb::Matrix::zeros(m, m);

      memset(M.data, 0, sizeof(double)*m*m);

      // Unpack symmetric.
      for (int i=0; i<(int)indices.size(); i++) {
	M(indices[i].first, indices[i].second) = data[i]; 	
	M(indices[i].second, indices[i].first) = data[i]; 
      }
    }
  }
};



/**
 */
class VariationalBasisProblemUV {
public:

  constexpr static double omega = 1.25;
  mg::restriction_image restrict;

  static Vec2d defaultBeta;
  static double defaultAlpha;

  Vec2d beta;  //!<Weight of smoothness term.
  std::vector<double>   S_alpha; // !<Weight of data terms
  std::vector<nacb::Imagef> S_I; // !<Data term tensors
  std::vector<MatrixPacker> S_I_packers;

  nacb::Imagef D_0, d;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  std::vector<nacb::Imagef> weights; //!<Per-pixel weighting, one for each S_I

  double              hx, hy;
  
  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);

  std::vector<int> times;
  std::vector<double> perFrameBeta;
  FlowBasis3D::ptr basis;

  VariationalBasisProblemUV(int width, int height){ 
    d_phi_D = d_phi_1en2;
    d_phi_S  = d_phi_1en2;    
  }

  VariationalBasisProblemUV(const std::vector<DisplaceOffsetUV> & warps, // Mapping from UV-space to image (for each image)
			    const std::vector<std::vector<Imagef> > & images,
			    FlowBasis3D::ptr& basisInput,
			    int referenceTimeIndex,
			    BaseMeshAnimation& baseMeshAnimation,
			    const std::vector<std::vector<nacb::Imagef> > & indyWeights,
			    std::vector<std::pair<int, int> > & pairs,
			    const nacb::Imagef & _D_0,
			    const std::vector<int>& timesInput,
			    const Vec2d & _alpha = Vec2d(1, 1),
			    const Vec2d & _beta = VariationalBasisProblemUV::defaultBeta,
			    double _hx = 1.0, double _hy = 1.0) : 
    beta(_beta), D_0(_D_0), hx(_hx), hy(_hy) {
    nacb::Timer timer;
      
    times = timesInput;
    basis = basisInput;

    // For now initialize the per-frame regularization with beta.y
    perFrameBeta = std::vector<double>(times.size(), _beta.y);

    int width = D_0.w, height = D_0.h;
    int numFrames = images[0].size();

    assert(indyWeights.size() == warps.size());
    
    d_phi_D = d_phi_1en2;
    d_phi_S  = d_phi_1en2;

    int numBasis = basis->getNumDimensions();

    f_h = nacb::Imagef(width, height, numBasis);
    d = nacb::Imagef(width, height, numBasis);
    residuals = nacb::Imagef(width, height, numBasis);
    
    int totalTensorParams = 0;
    int usedTensorParams = 0;

    // Stereo terms, one for each time
    for(int i=0; i<(int)pairs.size(); i++){
      for (int j=0; j<numFrames; j++) {
	if (basis->isSparse())
	  S_I_packers.push_back(MatrixPacker(basis->getSparseInteractions(j)));
	else
	  S_I_packers.push_back(MatrixPacker(numBasis, true));

	S_I.push_back(nacb::Imagef(width, height, S_I_packers.back().getPackedSize()));
	S_I.back() = 0;
	S_alpha.push_back(_alpha.x / numFrames); // Divide by num frame times.

	weights.push_back(indyWeights[pairs[i].first][j] * 
			  indyWeights[pairs[i].second][j]);

	totalTensorParams += (numBasis + 1)*(numBasis + 2)/2;
	usedTensorParams += S_I_packers.back().getPackedSize();
      }
    }

    
    // Add all the flow terms.
    //int indyBase = S_I.size();
    if (_alpha.y > 0) {
      for(int i=0; i<(int)images.size(); i++){
	for (int j=0; j<int(numFrames) - 1; j++) {
	  if (basis->isSparse())
	    S_I_packers.push_back(MatrixPacker(basis->getSparseInteractions(j) + basis->getSparseInteractions(j + 1)));
	  else
	    S_I_packers.push_back(MatrixPacker(numBasis, true));
	  
	  S_I.push_back(nacb::Imagef(width, height, S_I_packers.back().getPackedSize()));
	  S_I.back() = 0;
	  S_alpha.push_back(_alpha.y / numFrames);
	  
	  weights.push_back(indyWeights[i][j] * indyWeights[i][j+1]);
	  
	  totalTensorParams += (numBasis + 1)*(numBasis + 2)/2;
	  usedTensorParams += S_I_packers.back().getPackedSize();
	}
      }
    }
    printf("usedTensorParams: %d  (dense: %d)\n", usedTensorParams, totalTensorParams);


    // A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    d = 0;
    f_h = 0;
    
    // Add the stereo terms.
    int constraintIndex = 0;
    for(int pi=0; pi<(int)pairs.size(); pi++){
      for (int ft=0; ft<numFrames; ft++, constraintIndex++) {
	baseMeshAnimation.setTime(times[ft]);

	int index1 = pairs[pi].first;
	int index2 = pairs[pi].second;

	const DisplaceOffsetUV & warp1 = warps[index1];
	const DisplaceOffsetUV & warp2 = warps[index2];
	const nacb::Imagef & im1 = images[index1][ft];
	const nacb::Imagef & im2 = images[index2][ft];

	const int nchans = std::min(3, im2.nchannels);
	
	nacb::Imagef Ix1, Iy1;
	nacb::Imagef Ix2, Iy2;

	image_gradient_all(im1, Ix1, Iy1);
	image_gradient_all(im2, Ix2, Iy2);
	
	Ix1 *= (1.0/hx);
	Iy1 *= (1.0/hy);
	
	Ix2 *= (1.0/hx);
	Iy2 *= (1.0/hy);

#ifdef USE_OMP
#pragma omp parallel for
#endif
	for(int y=0; y<height; y++){
	  for(int x=0; x<width; x++){
	    double displace = 0;
	    Vec3f offs = (*basis)(displace, x, y, D_0, ft);
	  
	    nacb::Vec2d co1_o_deriv[3];
	    nacb::Vec2d co2_o_deriv[3];
	  
	    nacb::Vec2f co1 = warp1(x, y, displace, offs);
	    nacb::Vec2f co2 = warp2(x, y, displace, offs);
	    nacb::Vec2d co1_deriv;
	    nacb::Vec2d co2_deriv;
	  
	    co1_deriv = warp1.derivatives(x, y, displace, offs, co1_o_deriv);
	    co2_deriv = warp2.derivatives(x, y, displace, offs, co2_o_deriv);

	    // Take the derivative of the basis now
	    nacb::Matrix basisDerivMatrix = basis->derivative(x, y, D_0, ft);

	    for(int chan = 0; chan<nchans; chan++) {
	      Vec2d i1chan(Ix1.bilinear(co1.x, co1.y, chan),
			   Iy1.bilinear(co1.x, co1.y, chan));
	      
	      Vec2d i2chan(Ix2.bilinear(co2.x, co2.y, chan),
			   Iy2.bilinear(co2.x, co2.y, chan));
	      
	      double Iz = (im2.bilinear(co2.x, co2.y, chan) - 
			   im1.bilinear(co1.x, co1.y, chan));
	      
	      Vec4d i1chan_warp(i1chan.dot(co1_deriv),
				i1chan.dot(co1_o_deriv[0]),
				i1chan.dot(co1_o_deriv[1]),
				i1chan.dot(co1_o_deriv[2]));
	      
	      Vec4d i2chan_warp(i2chan.dot(co2_deriv),
				i2chan.dot(co2_o_deriv[0]),
				i2chan.dot(co2_o_deriv[1]),
				i2chan.dot(co2_o_deriv[2]));
	      // Take dot product of basisDeriv with image derivatives.
	      nacb::Matrix II(numBasis + 1, 1);

	      for (int i=0; i<numBasis; ++i) {
		nacb::Vec4d basisDeriv(basisDerivMatrix(0, i), basisDerivMatrix(1, i), 
				       basisDerivMatrix(2, i), basisDerivMatrix(3, i));
		II[i] = (dot4(i2chan_warp, basisDeriv) - dot4(i1chan_warp, basisDeriv));
	      }
	      II[numBasis] = Iz;
	      nacb::Matrix p = S_I_packers[constraintIndex].pack(II*II.transpose());
	      
	      for (int k=0; k<p.m*p.n; k++)
		S_I[constraintIndex](x, y, k) += p[k];
	    }
	  }
	}
      }
    }


    if (_alpha.y > 0) {

      for(int i=0; i<(int)images.size(); i++) {
	const DisplaceOffsetUV & warp = warps[i];
	const std::vector<nacb::Imagef> & im = images[i];
	const int nchans = std::min(3, im[0].nchannels);
	
      
	std::vector<nacb::Matrix> referenceVectors(nchans*width*height, nacb::Matrix(1, 1));

	// Compute the components from the reference frame (this is used in all others).
	{
	  nacb::Imagef Ix, Iy;
	  image_gradient_all(images[i][referenceTimeIndex], Ix, Iy);
	  Ix *= (1.0/hx);
	  Iy *= (1.0/hy);

	  baseMeshAnimation.setTime(referenceTimeIndex);

	  for(int y=0; y<height; y++){
	    for(int x=0; x<width; x++){
	      double displace = 0;
	      Vec3f offs = (*basis)(displace, x, y, D_0, referenceTimeIndex);
	      nacb::Vec2d co_o_deriv[3];
	      nacb::Vec2f co = warp(x, y, displace, offs);

	      nacb::Vec2d co_deriv;
	      co_deriv = warp.derivatives(x, y, displace, offs, co_o_deriv);

	      // Take the derivative of the basis now
	      nacb::Matrix basisDerivMatrix = basis->derivative(x, y, D_0, referenceTimeIndex);
	    
	      for(int chan = 0; chan<nchans; chan++){
		Vec2d i0chan(Ix.bilinear(co.x, co.y, chan),
			     Iy.bilinear(co.x, co.y, chan));

		Vec4d i0chan_warp(i0chan.dot(co_deriv),
				  i0chan.dot(co_o_deriv[0]),
				  i0chan.dot(co_o_deriv[1]),
				  i0chan.dot(co_o_deriv[2]));

		// Take dot product of basisDeriv with image derivatives.
		nacb::Matrix II(numBasis + 1, 1);

		for (int i=0; i<numBasis; ++i) {
		  nacb::Vec4d basisDeriv(basisDerivMatrix(0, i), basisDerivMatrix(1, i), 
					 basisDerivMatrix(2, i), basisDerivMatrix(3, i));
		  II[i] = dot4(i0chan_warp, basisDeriv);
		}
		II[numBasis] = im[referenceTimeIndex].bilinear(co.x, co.y, chan);

		referenceVectors[nchans*(y*width + x) + chan] = II;
	      }
	    }
	  }	
	}

	for (int ft=0; ft<numFrames; ft++, constraintIndex++) {
	  if (ft == referenceTimeIndex) {
	    constraintIndex--;
	    continue;
	  }
	  baseMeshAnimation.setTime(ft);

	  nacb::Imagef Ix, Iy;
	  image_gradient_all(images[i][ft], Ix, Iy);
	  Ix *= (1.0/hx);
	  Iy *= (1.0/hy);

	
#ifdef USE_OMP
#pragma omp parallel for
#endif
	  for(int y=0; y<height; y++){
	    for(int x=0; x<width; x++){
	      double displace = 0;
	      Vec3f offs = (*basis)(displace, x, y, D_0, ft);
	      nacb::Vec2d co_o_deriv[3];
	      nacb::Vec2f co = warp(x, y, displace, offs);

	      nacb::Vec2d co_deriv;
	      co_deriv = warp.derivatives(x, y, displace, offs, co_o_deriv);

	      // Take the derivative of the basis now
	      nacb::Matrix basisDerivMatrix = basis->derivative(x, y, D_0, ft);
	    
	      for(int chan = 0; chan<nchans; chan++){
		Vec2d i2chan(Ix.bilinear(co.x, co.y, chan),
			     Iy.bilinear(co.x, co.y, chan));

		Vec4d i2chan_warp(i2chan.dot(co_deriv),
				  i2chan.dot(co_o_deriv[0]),
				  i2chan.dot(co_o_deriv[1]),
				  i2chan.dot(co_o_deriv[2]));

		// Take dot product of basisDeriv with image derivatives.
		nacb::Matrix II(numBasis + 1, 1);

		for (int i=0; i<numBasis; ++i) {
		  nacb::Vec4d basisDeriv(basisDerivMatrix(0, i), basisDerivMatrix(1, i), 
					 basisDerivMatrix(2, i), basisDerivMatrix(3, i));
		  II[i] = (dot4(i2chan_warp, basisDeriv) - referenceVectors[nchans*(y*width + x) + chan][i]);
		}
		II[numBasis] = im[ft].bilinear(co.x, co.y, chan)
		  - referenceVectors[nchans*(y*width + x) + chan][numBasis];

		nacb::Matrix p = S_I_packers[constraintIndex].pack(II*II.transpose());
	      
		for (int k=0; k<p.m*p.n; k++)
		  S_I[constraintIndex](x, y, k) += p[k];
	      }
	    }
	  }
	}
      }
    }

    std::cout << "Num tensors:" << S_I.size() << std::endl;
    std::cout << "Setting up problem:" << double(timer) << std::endl;

    // DisplaceOffsetBasis::normalizeWeights(weights);
  }
 
  VariationalBasisProblemUV  restriction();
  
  nacb::Imagef prolong(const nacb::Imagef & img){
    return restrict.prolong(img);
  }
  
  robust_stats iteration(nacb::Imagef & residuals, bool update = true);
   
  void smooth(){
    iteration(residuals, true);
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
  static double imageDataEnergy(const nacb::Imagef & D,
				const std::vector<std::vector<nacb::Imagef> > & images,
				const std::vector<DisplaceOffsetUV> & warps,
				const std::vector<std::vector<nacb::Imagef> > & weights,
				const std::vector<std::pair<int, int> > & pairs,
				const std::vector<int>& times,
				int referenceTime,
				const FlowBasis3D::ptr & basis,
				BaseMeshAnimation& baseMeshAnimation,
				double (* phi_D)(double ) = phi_1en2){
    nacb::Timer timer;

    int width = D.width, height = D.height;
    double energy_D = 0;
    int numFrames = images[0].size();

    baseMeshAnimation.setTime(times[referenceTime]);
    std::vector<nacb::Imagef> referenceColors(images.size());

    for (int i=0; i<(int)images.size(); ++i) {
      nacb::Imagef colors(width, height, 4);
      referenceColors[i] = colors;

      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  double displace = 0;
	  Vec3f offs = (*basis)(displace, x, y, D, referenceTime);
	  
	  nacb::Vec2f co1 = warps[i](x, y, displace, offs);
	  images[i][referenceTime].bilinear(co1.x, co1.y, &(colors(x, y, 0)));
	}
      }
    }

    for (int ft=0; ft<(int)images[0].size(); ft++) {
      if (ft == referenceTime) continue;
      
      baseMeshAnimation.setTime(times[ft]);

      for (int i=0; i<(int)images.size(); i++) {
	for(int y=0; y<height; y++){
	  for(int x=0; x<width; x++){
	    double displace = 0;
	    Vec3f offs = (*basis)(displace, x, y, D, ft);
	    nacb::Vec2f co = warps[i](x, y, displace, offs);
	    nacb::Vec3f color;
	    
	    images[i][ft].bilinear(co.x, co.y, color.data, 3);
	    nacb::Vec3f diff = color - nacb::Vec3f(referenceColors[i](x, y, 0),
						   referenceColors[i](x, y, 1),
						   referenceColors[i](x, y, 2));

	    energy_D += weights[i][ft](x, y)*phi_D(diff.dot(diff));
	  }
	}
      }
    }


    for (int i=0; i<(int)pairs.size(); i++) {
      int i1 = pairs[i].first, i2 = pairs[i].second;
      
      for (int ft=0; ft < int(numFrames); ft++) {
	baseMeshAnimation.setTime(times[ft]);
		
	for(int y=0; y<height; y++){
	  for(int x=0; x<width; x++){	
	    double displace = 0;
	    Vec3f offs = (*basis)(displace, x, y, D, ft);

	    nacb::Vec2f co1 = warps[i1](x, y, displace, offs);
	    nacb::Vec2f co2 = warps[i2](x, y, displace, offs);
	    nacb::Vec3f color1, color2;
	    
	    images[i1][ft].bilinear(co1.x, co1.y, color1.data, 3);
	    images[i2][ft].bilinear(co2.x, co2.y, color2.data, 3);
	    
	    nacb::Vec3f diff = color2 - color1;
	    energy_D += weights[i1][ft](x, y)*weights[i2][ft](x, y)*phi_D(diff.dot(diff));
	  }
	}
      }
    }
    std::cout << "Computing energy took:" << double(timer) << std::endl;
    return energy_D;
  }

  static double smoothnessEnergy(const nacb::Imagef & D,
				 double (* phi_S)(double ) = phi_1en2){
    double energy_S = 0;
    int width = D.width;
    int height = D.height;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	int xup = std::min(x + 1, width - 1);
	int yup = std::min(y + 1, height - 1);
	int xdn = std::max(x - 1, 0);
	int ydn = std::max(y - 1, 0);
	
	nacb::Vec2d grad((D(xup, y) - D(xdn, y))/2.0,
			 ((D(x, yup) - D(x, ydn))/2.0));

	energy_S += phi_S(grad.dot(grad));
      }
    }
    return energy_S;
  }

  static double smoothnessEnergyFlow(const nacb::Imagef & D,
				     double (* phi_S)(double ) = phi_1en2){
    double energy_S = 0;
    int width = D.width;
    int height = D.height;

    nacb::Imagef gs[3] = {D.gradient(1), D.gradient(2), D.gradient(3)};

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	nacb::Vec2d g0(gs[0](x, y, 0), gs[0](x, y, 1));
	nacb::Vec2d g1(gs[1](x, y, 0), gs[1](x, y, 1));
	nacb::Vec2d g2(gs[2](x, y, 0), gs[2](x, y, 1));
	
	energy_S += phi_S(g0.dot(g0) + g1.dot(g1) + g2.dot(g2));
      }
    }
    return energy_S;
  }

};



#endif
