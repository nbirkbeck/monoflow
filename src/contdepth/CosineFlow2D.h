/** \brief I have been thinking about using a cosine basis for modeling the 
           motion of a surface in 3D over time.

  The idea is that in 3D these constraints will allow reconstruction from
  intra-camera observations.  This 2D example is just a test.  It works
  with several different basis types.

  The objective is currently:
  \sum_{i=1}^N\int phi(|I_i(x(i)) - I_0(x(0))|^2) + \sum \int phi(|\nabla b_j|^2)
   
  x(i) = \sum b_j B_j(i)

  The 2D basis elements either affect the x or y direction (independently).
  There are linear basis vectors (always 2 basis elements) or the cosine basis.

  The cosine basis also uses a constraint to keep the flow at time 0 consistent:
  phi(|x(0)|^2) + phi(|y(0)|^2)

  This actually causes the cosine basis to be worse in simple 2 frame flow.  
  I tested these ideas with some of the data sets from middlebury and they seemed
  to work fairly well.

  ISSUES: -the normalization of the basis functions is currently done 
          for each basis image independently.  This could easily be
          combined into say phi(|\nabla b_1|^2 + |\nabla b_2|^2)
          -parameter settings.
          -I noticed that the optimizer had trouble when using phi_1en6;
           it seemed to bounce around.  Using phi_1en2 converges better. 
           I don't think this is an issue with the code in "iteration",
           although there could be a bug.
          -It might be smart to use some sort of initialization for the cosine
           basis (e.g., compute two frame flow and project), or at least keep
           increasing the weihgt for later frames (or introduce later frames
           slowly).

  Neil Birkbeck
*/
#ifndef COSINE_FLOW_H
#define COSINE_FLOW_H

#include "autorecon/stereo/flow_variational.h"
#include "autorecon/stereo/multigrid.h"

#include "Math.h"
#include <nmisc/timer.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>

#include <nimage/image.h>
#include <vector>
#include <stdexcept>

/** \brief Evaluate the basis densely to produce a flow map.
 */
template <class BasisType>
nacb::Imagef basisToFlow(BasisType& basis, nacb::Imagef& D, int t) {
  nacb::Imagef flow(D.w, D.h, 2);

  for (int y=0; y<D.h; y++) {
    for (int x=0; x<D.w; x++) {
      nacb::Vec2f flowed = basis(x, y, D, t);
      flow(x, y, 0) = flowed.x - x;
      flow(x, y, 1) = flowed.y - y;
    }
  }

  return flow;
};


/** \brief This is the standard optic-flow basis function.
 */
class StandardFlowBasis {
 public:
  int getNumDimensions() const { return 2; }
  
  nacb::Vec2f operator()(float x, float y, const nacb::Imagef& D, int t) {
    return nacb::Vec2f(x + D(x, y, 0)*t, y + D(x, y, 1)*t);
  }
  
  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    return nacb::Matrix::eye(2, 2)*t;
  }

  bool needsConstraintAtInitialTime() { return false; }
};


/** \brief A test of just translation in the x-direction.
 */
class FlowBasisX {
 public:
  int getNumDimensions() const { return 1; }
  
  nacb::Vec2f operator()(float x, float y, const nacb::Imagef& D, int t) {
    return nacb::Vec2f(x + D(x, y, 0)*t, y);
  }
  
  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    return nacb::Matrix::eye(2, 1)*t;
  }

  bool needsConstraintAtInitialTime() { return false; }
};


/** \brief A cosine basis for the motion.
 */
class CosineBasis {
 public:
  CosineBasis(int N):m_N(2*N) { }

  int getNumDimensions() {return m_N; }

  nacb::Vec2f operator()(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Vec2f dx(x, y);
    for (int i=0; i<D.nchannels; i++) {
      int basisIndex = int(i/2);
      double value = cos(M_PI/m_N * (0.5 + t)*basisIndex) * D(x, y, i);
      dx.data[i % 2] += value;
    }
    return dx;
  }
  
  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Matrix deriv = nacb::Matrix::zeros(2, m_N);

    for (int i=0; i<D.nchannels; i++) {
      int basisIndex = int(i/2);
      double value = cos(M_PI/m_N * (0.5 + t)*basisIndex);
      deriv(i % 2, i) = value;
    }
    return deriv;
  }

  bool needsConstraintAtInitialTime() { return true; }
 protected:
  int m_N;
};


/** \brief The cosine flow takes a set of basis functions and tries to do optic
           flow over several frames.
	   

   The cosine basis needed constraints to keep first frame offsets zero.  This
   caused the results to be slightly worse than when using the linear basis.
 */
template <class WarpType>
class CosineFlow2D {
public:

  constexpr static double omega = 1.25;
  mg::restriction_image restrict;

  std::vector<double>   beta; // !<Weight of smooth terms.
  std::vector<double>   S_alpha; // !<Weight of data terms
  std::vector<nacb::Imagef> S_I; // !<Data term tensors

  nacb::Imagef D_0, d;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  double              hx, hy;
  
  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);


  CosineFlow2D(int width, int height){ 
    d_phi_D = d_phi_1en2;
    d_phi_S  = d_phi_1en2;    
  }

  CosineFlow2D(const std::vector<nacb::Imagef> & images, 
	       WarpType& warp,
	       const nacb::Imagef & _D_0,
	       double _alpha = 1.0,
	       double _beta = 6.0/255.0,
   	       double _hx = 1.0, double _hy = 1.0) : 
    beta(warp.getNumDimensions(), _beta), D_0(_D_0), hx(_hx), hy(_hy) {

    int width = D_0.w, height = D_0.h;
    
    d_phi_D = d_phi_1en2;
    d_phi_S  = d_phi_1en2;

    f_h = nacb::Imagef(width, height, warp.getNumDimensions());
    d = nacb::Imagef(width, height, warp.getNumDimensions());
    residuals = nacb::Imagef(width, height, warp.getNumDimensions());

    for(int i=1; i<(int)images.size(); i++){
      int n = warp.getNumDimensions() + 1;
      S_I.push_back(nacb::Imagef(width, height, n*(n + 1)/2));
      S_alpha.push_back(_alpha);
      
      S_I[i - 1] = 0;
    }

    // A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    d = 0;
    f_h = 0;
    

    for(int i=1; i<(int)images.size(); i++){
      int index1 = 0;
      int index2 = i;

      const nacb::Imagef & im1 = images[index1];
      const nacb::Imagef & im2 = images[index2];
      
      const int nchans = std::min(3, im2.nchannels);
      
      // Assume that the images are not warped.
      nacb::Imagef im1_grad[3] = {im1.gradient(0),
				  im1.gradient(1),
				  im1.gradient(2)};

      nacb::Imagef im2_grad[3] = {im2.gradient(0),
				  im2.gradient(1),
				  im2.gradient(2)};

      nacb::Imagef Ix1[3] = {im1_grad[0].getChannel(0)*(1.0/hx),
			     im1_grad[1].getChannel(0)*(1.0/hx),
			     im1_grad[2].getChannel(0)*(1.0/hx)};
			    
      nacb::Imagef Iy1[3] = {im1_grad[0].getChannel(1)*(1.0/hy),
			     im1_grad[1].getChannel(1)*(1.0/hy),
			     im1_grad[2].getChannel(1)*(1.0/hy)};

      nacb::Imagef Ix2[3] = {im2_grad[0].getChannel(0)*(1.0/hx),
			     im2_grad[1].getChannel(0)*(1.0/hx),
			     im2_grad[2].getChannel(0)*(1.0/hx)};
			    
      nacb::Imagef Iy2[3] = {im2_grad[0].getChannel(1)*(1.0/hy),
			     im2_grad[1].getChannel(1)*(1.0/hy),
			     im2_grad[2].getChannel(1)*(1.0/hy)};

      nacb::Matrix Ix1_sample(1, 2);
      nacb::Matrix Ix2_sample(1, 2);

      nacb::Matrix Iy1_sample(1, 2);
      nacb::Matrix Iy2_sample(1, 2);

      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){

	  // Already displaced coordiantes.
	  nacb::Vec2f co1 = warp(x, y, _D_0, index1);
	  nacb::Vec2f co2 = warp(x, y, _D_0, index2);
	  nacb::Matrix co1_deriv = warp.derivative(x, y, _D_0, index1);
	  nacb::Matrix co2_deriv = warp.derivative(x, y, _D_0, index2);
	
	  for(int chan = 0; chan<nchans; chan++){
	    Ix1_sample[0] = Ix1[chan].bilinear(co1.x, co1.y);
	    Ix1_sample[1] = Iy1[chan].bilinear(co1.x, co1.y);

	    Ix2_sample[0] = Ix2[chan].bilinear(co2.x, co2.y);
	    Ix2_sample[1] = Iy2[chan].bilinear(co2.x, co2.y);

	    nacb::Matrix der1 = Ix1_sample * co1_deriv;
	    nacb::Matrix der2 = Ix2_sample * co2_deriv;
	      
	    double Iz = im2.bilinear(co2.x, co2.y, chan) - im1.bilinear(co1.x, co1.y, chan);
	    nacb::Matrix der = der2 - der1;

	    // n = warp.getNumDimensions() + 1)
	    // S is of size n * (n + 1) / 2
	    nacb::Matrix derIz(1, der.n + 1);
	    derIz.setSubmatrix(0, 0, der);
	    derIz[der.n] = Iz;

	    nacb::Matrix p = packSymmetricMatrix(derIz.transpose()*derIz);
	    if (S_I[i - 1].nchannels != p.m*p.n) {
	      std::cerr << "Invalid matrix packing.h:" << S_I[i - 1].nchannels << " " << p.m << " " << p.n << std::endl;
	      throw std::runtime_error("invalid matrix packing.");
	    }

	    for (int k=0; k<p.m*p.n; k++)
	      S_I[i - 1](x, y, k) += p[k];
	  }
	}
      }
    }

    // Add a constraint such that x^2 = 0 at time = 0 (and for y)
    if (warp.needsConstraintAtInitialTime()) {
      int n = warp.getNumDimensions() + 1;
      
      int coordConstraints[2] = {(int)S_I.size(), (int)S_I.size() + 1};
      for (int i=0; i<2; i++) {
	S_I.push_back(nacb::Imagef(width, height, n*(n + 1)/2));
	S_alpha.push_back(10);
	S_I.back() = 0;
      }
      
      for(int y=0; y<height; y++){
	for(int x=0; x<width; x++){
	  nacb::Vec2f co1 = warp(x, y, _D_0, 0);
	  nacb::Matrix co1_deriv = warp.derivative(x, y, _D_0, 0);
	  co1.x -= x;
	  co1.y -= y;
	  
	  // Add a constraint for x and y
	  for (int i=0; i<2; i++) {
	    nacb::Matrix der = co1_deriv.getRow(i);
	    
	    nacb::Matrix derIz(1, der.n + 1);
	    derIz.setSubmatrix(0, 0, der);
	    derIz[der.n] = co1.data[i];
	    
	    nacb::Matrix p = packSymmetricMatrix(derIz.transpose()*derIz);
	    
	    for (int k=0; k<p.m*p.n; k++)
	      S_I[coordConstraints[i]](x, y, k) += p[k];
	  }
	}
      }
    }
  }
 
  CosineFlow2D  restriction();
  
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
				WarpType& warp,
				const std::vector<nacb::Imagef> & images,
				double (* phi_D)(double ) = phi_1en6){
    int width = D.width, height = D.height;
    double energy_D = 0;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){	
	for (int i=1; i<(int)images.size(); i++) {
	  int i1 = 0, i2 = i;
	  
	  nacb::Vec2f co1 = warp(x, y, D, i1);
	  nacb::Vec2f co2 = warp(x, y, D, i2);
	  nacb::Vec3f color1, color2;
	  
	  images[i1].bilinear(co1.x, co1.y, color1.data, 3);
	  images[i2].bilinear(co2.x, co2.y, color2.data, 3);
	  
	  nacb::Vec3f diff = color2 - color1;
	  energy_D += phi_D(diff.dot(diff));
	}
      }
    }
    return energy_D;
  }

  static double smoothnessEnergy(const nacb::Imagef & D,
				 double (* phi_S)(double ) = phi_1en6){
    double energy_S = 0;
    int width = D.width;
    int height = D.height;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
  	int neigh[4][2] = {{-1,0}, {1,0}, {0,-1}, {0, 1}};
	double hs[4] = {1, 1, 1, 1};
	
	int xup = std::min(x + 1, width - 1);
	int yup = std::min(y + 1, height - 1);
	int xdn = std::max(x - 1, 0);
	int ydn = std::max(y - 1, 0);
	
	for (int k=0; k<D.nchannels; k++) {
	  nacb::Vec2d grad((D(xup, y, k) - D(xdn, y, k))/2.0,
			   ((D(x, yup, k) - D(x, ydn, k))/2.0));
	  
	  energy_S += phi_S(grad.dot(grad));
	}
      }
    }
    return energy_S;
  }
};



#endif
