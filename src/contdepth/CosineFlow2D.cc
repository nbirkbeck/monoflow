#include "CosineFlow2D.h"
#include <nmisc/commandline.h>
#include "autorecon/stereo/flow.h"

using namespace std;
using namespace nacb;
using namespace mg;



template <class WarpType>
CosineFlow2D<WarpType> CosineFlow2D<WarpType>::restriction(){
  int neww = (int)ceil(get_width()/2), newh = (int)ceil(get_height()/2);
  
  restrict = restriction_image(get_width(), get_height(), neww, newh);
    
  CosineFlow2D sub(neww, newh);
  sub.S_I.clear();

  for(int i=0; i<(int)S_I.size(); i++){
    sub.S_I.push_back(restrict.apply(S_I[i]));
  }
  sub.S_alpha = S_alpha;

  sub.D_0  = restrict.apply(D_0);
  sub.d    = restrict.apply(d);

  sub.hx =  hx * double(get_width())/double(neww);
  sub.hy =  hy * double(get_height())/double(newh);

  sub.beta = beta;
  sub.residuals = nacb::Imagef(neww, newh, D_0.nchannels);


  nacb::Imagef residuals_H;
  sub.f_h = nacb::Imagef(neww, newh, D_0.nchannels);
  sub.f_h = 0;
  sub.iteration(residuals_H, false);
    
  //residuals_H = 0 -AH*uv_H  => AH*uv_H = -residuals_H
  sub.f_h = restrict.apply(residuals) - residuals_H;// = r_H +  AH*uv_H;
  return sub;
}


template <class WarpType>
robust_stats CosineFlow2D<WarpType>::iteration(nacb::Imagef & residuals, //Input output
					       bool update){
  nacb::StopWatch setup;

  nacb::Imagef d_old = d.copy();
  nacb::Imagef D = D_0 + d;
  std::vector<nacb::Imagef> D_phi_S(D_0.nchannels);

  // bnorm: Each channel of the basis is regularized independently.
  // bnorm: If you just change D_phi_s here to merge what terms you want.
  // bnorm: Also change classInd below.
  for (int i=0; i<D_0.nchannels; ++i) {
    nacb::Imagef D_grad = D.gradient(i);
    D_phi_S[i] = get_D_phi_S(D_grad, hx, hy, d_phi_S);
  }

  //printf("setup took %lf\n", double(setup));
  double res = 0.0;
  double len = 0.0;
  double moved = 0.0;

  int width = S_I[0].w, height = S_I[0].h;
  int nchans = S_alpha.size();

  if(residuals.w != width || residuals.h != height || residuals.nchannels != 1){
    residuals = nacb::Imagef(width, height, D_0.nchannels);
    residuals = 0;
  }

  nacb::Matrix d_xy_vec(D_0.nchannels + 1, 1);

  // Store  matrices outside so they don't need to be reallocated.
  nacb::Matrix M(D_0.nchannels, D_0.nchannels);
  nacb::Matrix rhs(D_0.nchannels, 1);
  nacb::Matrix Msum(D_0.nchannels, 1);
  nacb::Matrix d_mat(D_0.nchannels, 1);
  nacb::Matrix upd(D_0.nchannels, 1);
  nacb::Matrix diff(D_0.nchannels, 1);
  nacb::Matrix d_S_xy(D_0.nchannels, 1);
  nacb::Matrix S[nchans];
  nacb::Matrix quad_form_temp;
  nacb::Matrix S_dd_xy_temp;
  
  nacb::Matrix dd_xy(D_0.nchannels + 1, 1);
  nacb::Matrix d_T_xy(D_0.nchannels + 1, 1);
  nacb::Matrix res_vec(D_0.nchannels, 1);


  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){	
      // Image data term
      double d_phi_D_xy[nchans];

      // Current solution
      for (int i=0; i<D_0.nchannels; ++i)
	d_xy_vec[i] = d_old(x, y, i);
      d_xy_vec[D_0.nchannels] = 1;

      for(int si=0; si<nchans; si++){
	unpackSymmetricMatrix<float>(&(S_I[si](x, y, 0)), S_I[si].nchannels, S[si]);
	//Vec2d(S_I[si](x, y, 0), S_I[si](x, y, 1));//unpack_sym_mat2x2(S_I[si], x, y);
	//Vec2d Ssi2 = Vec2d(S_I[si](x, y, 1), S_I[si](x, y, 2));
	//Vec2d Ssi_dx(S[si].dot(d_xy_vec), Ssi2.dot(d_xy_vec));
	d_phi_D_xy[si] = d_phi_D(quadratic_form_eval(S[si], d_xy_vec, quad_form_temp));
      }
	
      int neigh[4][2] = {{-1,0}, {1,0}, {0,-1}, {0, 1}};
      double hs[4] = {hx, hx, hy, hy};


      // !<per-pixel weighting function (on the data terms only).
      double wt_xy = 1.0;// weight(x, y);

      M.setAll(0);

      for(int k=0; k<rhs.m; k++){
	rhs[k] = f_h(x, y, k); 
      }

      for(int si=0; si<(int)S_I.size(); si++){
	for(int i=0; i<D_0.nchannels; i++){
	  for(int j=0; j<D_0.nchannels; j++){
	    M(i, j) += S[si](i, j)*(d_phi_D_xy[si]*S_alpha[si]*wt_xy);
	  }
	  rhs[i] -=  S[si](i, D_0.nchannels)*(d_phi_D_xy[si]*S_alpha[si]*wt_xy);
	}
	// + alpha*d_phi_D2_xy*T(0, 0);	
      }

      Msum.setAll(0);
      d_S_xy.setAll(0);

      for(int ni=0; ni<4; ni++){	
	int xn = x + neigh[ni][0];
	int yn = y + neigh[ni][1];
	  
	// Reflect the boundaries.
	if(xn<0)xn = 0;
	if(yn<0)yn = 0;
	if(xn>=width)xn = width-1;
	if(yn>=height)yn = height-1;
	  
	float * dd_xnyn = &d(xn, yn, 0);
	
	// Use most up to date solution.
	for(int j=0; j<D_0.nchannels; j++) {
	  int class_ind = j; // bnorm: Regularize each one separately.
	  double d_phi_S_xy = D_phi_S[class_ind](x, y);
	  double d_phi_S_xnyn = D_phi_S[class_ind](xn, yn);
	  
	  double d_phi_S_avg = (d_phi_S_xy + d_phi_S_xnyn)/2.0;
	  Msum[j] += beta[class_ind]*d_phi_S_avg/(hs[ni]*hs[ni]);

	  // What would be used in the actual equations.
	  d_S_xy[j] += -beta[class_ind] * ((d_phi_S_xy + d_phi_S_xnyn)/2.0 * 
					   (D(xn, yn, j) - D(x, y, j))/(hs[ni]*hs[ni]));
	  
	  
	  rhs[j] += beta[class_ind]*d_phi_S_avg*(D_0(xn, yn, j) + dd_xnyn[j] - D_0(x, y, j))/(hs[ni]*hs[ni]);
	}
      }
      
      for(int j=0; j<D_0.nchannels; j++) {
	M(j, j) += Msum[j];

	d_mat[j] = d(x, y, j);
      }
      
      if(isNormal(M) && isNormal(rhs)) {
	solve_symmetric_nxn(M, rhs, upd);

	if (1e4 < upd.dot(upd)) {
	  M.printMatlab("M");
	  rhs.printMatlab("rhs");
	  upd.printMatlab("upd");
	  copy_matrix_no_alloc(upd, d_mat);
	  printf("\n\n");
	}
      }
      else
	copy_matrix_no_alloc(upd, d_mat);

      if(!isNormal(M))
	copy_matrix_no_alloc(upd, d_mat);

      diff_matrix_no_alloc(diff, upd, d_mat);
            
      if(update){
	for(int j=0; j<D_0.nchannels; j++){
	  // Update solution, omega allows for SOR-type solution
	  d(x, y, j) = d(x, y, j)*(1.0-omega) + omega*upd[j];
	}
      }
      moved += diff.dot(diff);
      
      // These quantities would be used in the actual equations, just checking residual here.
      dd_xy[D_0.nchannels] = 1;
      for (int k=0; k<D_0.nchannels; ++k)
	dd_xy[k] = d(x, y, k);

      d_T_xy.setAll(0);

      for(int si=0; si<nchans; si++) {
	// S_dd_xy_temp = S[si]*dd_xy
	matrix_vec_mult(S[si], dd_xy, S_dd_xy_temp);
	d_T_xy += (S_dd_xy_temp)*(S_alpha[si]*wt_xy*d_phi_D_xy[si]);
      }
      
      for(int si=0; si<D_0.nchannels; si++) {
	res_vec[si] = (f_h(x, y, si) - (d_T_xy[si] + d_S_xy[si]));
      }
      // Should work also, except M has been modified by the solver...
      //double r = rhs - (M*d(x, y));
      // printf("%f %f\n", r, res_0);
      res += res_vec.dot(res_vec);
      len += upd.dot(upd);
      
      // These should probably use the updated coefficients, e.g., be recomputed after the loop,
      // when all the du, dv displacements have been updated.
      for (int si=0; si<D_0.nchannels; si++)
	residuals(x, y, si) = res_vec[si];
    }
  }
  //printf("%f  %f  moved: %f\n", res, len, moved);
  double erel = res/len;  
  return robust_stats(moved, res, len, erel);
}



template <class BasisType> 
void run_basis_flow(BasisType& basis, 
		    std::vector<nacb::Imagef> & images,
		    int downsample,
		    int maxIterations,
		    double alpha = 1.0)
{
  nacb::Imagef D(images[0].w, images[0].h, basis.getNumDimensions());
  D = 0;

  double imageEnergy = CosineFlow2D<BasisType>::imageDataEnergy(D, basis, images);
  std::cout << "image energy:" << imageEnergy << std::endl;

  for (int ds=downsample; ds>=0; ds--) {
    std::vector<nacb::Imagef> imagesDownsample = images;

    if (ds != 0) {
      imagesDownsample.clear();
      for (int i=0; i<(int)images.size(); ++i)
	imagesDownsample.push_back(images[i].resize(images[i].w / (1 << ds), 
						    images[i].h / (1 << ds)));
    }

    double factor = double(imagesDownsample[0].w)/D.w;
    D = D.resize(imagesDownsample[0].w, imagesDownsample[0].h) * factor;
    printf("%dx%d  (factor: %f, downsample: %d)\n", D.w, D.h, factor, ds);

    for (int its=0; its<maxIterations; its++) {
      CosineFlow2D<BasisType> csflow(imagesDownsample, basis, D, alpha/(1 << downsample));
      
      nacb::Imagef res;
      csflow.iteration(res, false).print();
      
      int nsolve = 1;
      mg::wcycle_solver solver(2, 2, 5, 5, true);
      //mg::basic_solver solver;nsolve = 100;
      
      solver.solve(csflow, nsolve);
      csflow.iteration(res, false).print();
      
      D = csflow.getCombinedResult();
      
      double imageEnergy = CosineFlow2D<BasisType>::imageDataEnergy(D, basis, imagesDownsample);
      std::cout << "image energy:" << imageEnergy << std::endl;
    }
  }

  // Get the min over all displacements.
  double minValue = 1e5, maxValue = -1e5;
  for (int i=0; i<(int)images.size(); i++) {
    double mn, mx;
    nacb::Imagef D0 = basisToFlow(basis, D, i);
    D0.getRange(mn, mx);
    printf("Flow %d range: %lf %lf\n", i, mn, mx);
    minValue = std::min(minValue, mn);
    maxValue = std::max(maxValue, mx);

  }

  // Output some images.
  for (int i=0; i<(int)images.size(); i++) {
    char name[1024];

    nacb::Imagef im = basisToFlow(basis, D, i);
    nacb::Vec2d avg(0, 0);
    for (int y=0; y<im.h; y++) {
      for (int x=0; x<im.w; x++) {
	avg.x += im(x, y, 0);
	avg.y += im(x, y, 1);
      }
    }
    printf("avg %d: %lf %lf (outputting to /tmp/dx* /tmp/dy* /tmp/warped-*)", i, avg.x/(im.w*im.h), avg.y/(im.w*im.h));
    
    snprintf(name, 1024, "/tmp/dx-%d.png", i);
    ((im.getChannel(0) - minValue)*(1.0/(maxValue - minValue))).save(name);

    snprintf(name, 1024, "/tmp/dy-%d.png", i);
    ((im.getChannel(1) - minValue)*(1.0/(maxValue - minValue))).save(name);

    nacb::Imagef warped = flow_warp_back(im, images[i]);
    snprintf(name, 1024, "/tmp/warped-%02d.png", i);
    warped.save(name);
  }
}


int main(int ac, char * av[]){
  std::string basisType = "cosine";
  int numBasis = 0;
  int downsample = 0;
  int maxIterations = 1;

  nacb::CommandLine cline("Run with a list of input images.");

  cline.registerOption("basis", "The basis type", &basisType);
  cline.registerOption("numBasis", "The number of basis elements", &numBasis);
  cline.registerOption("maxIterations", "The number of iterations", &maxIterations);
  cline.registerOption("downsample", "The number of times to downsample", &downsample);
  cline.parse(ac, av);

  std::cout << "Using basis:" << basisType << std::endl;
  std::vector<nacb::Imagef> images;

  for (int i=optind; i<ac; i++) {
    images.push_back(nacb::Imagef(av[i]));
  }

  if (basisType == "linear") {
    StandardFlowBasis basis;
    run_basis_flow(basis, images, downsample, maxIterations);
  }
  else {
    if (numBasis <= 0) numBasis = images.size();
    std::cout << "Num basis:" << numBasis << std::endl;
    CosineBasis basis(numBasis);
    run_basis_flow(basis, images, downsample, maxIterations);
  }

  return 0;
}
