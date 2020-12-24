#include "autorecon/stereo/flow_variational.h"
#include "VariationalBasisProblemUV.h"
#include "Math.h"

using namespace std;
using namespace nacb;
using namespace mg;


Vec2d VariationalBasisProblemUV::defaultBeta(0.5, 0.5);
double VariationalBasisProblemUV::defaultAlpha = 1.0;


nacb::Imagef get_D_phi_S(std::vector<nacb::Imagef>& D_grad, 
			 double hx, double hy,
			 double (*d_phi_S)(double),
			 const std::vector<int>& mapping,
			 int numMaps) {
  nacb::Imagef D_phi_S(D_grad[0].w, D_grad[0].h, numMaps);
  D_phi_S = 0;

  for(int y=0; y<D_phi_S.h; y++){
    for(int x=0; x<D_phi_S.w; x++){ 
      // Accumulate everything that is robustly regularized together.
      for (int k=0; k<(int)D_grad.size(); k++) {
	D_phi_S(x, y, mapping[k]) += (sqr(D_grad[k](x, y, 0)/hx) + sqr(D_grad[k](x, y, 1)/hy));
      }

      // Apply the function.
      for (int k=0; k<numMaps; k++) 
	D_phi_S(x, y, k) = d_phi_S(D_phi_S(x, y, k));
    }
  }
  return D_phi_S;
}


nacb::Imagef get_D_phi_S_perFrame(nacb::Imagef& D, 
				  double hx, double hy,
				  FlowBasis3D::ptr& basis,
				  std::vector<int>& times, 
				  double (*d_phi_S)(double)) {
  nacb::Imagef D_phi_S(D.w, D.h, times.size());
  nacb::Imagef offsets(D.w, D.h, 3);

  for (int t=0; t<(int)times.size(); t++) {
    for(int y=0; y<D_phi_S.h; y++){
      for(int x=0; x<D_phi_S.w; x++){ 
	double displace = 0;
	nacb::Vec3f offs = (*basis)(displace, x, y, D, times[t]);
	offsets(x, y, 0) = offs.x;
	offsets(x, y, 1) = offs.y;
	offsets(x, y, 2) = offs.z;
      }
    }
    
    nacb::Imagef G[3] = {offsets.gradient(0),
			 offsets.gradient(1),
			 offsets.gradient(2)};

    for(int y=0; y<D_phi_S.h; y++){
      for(int x=0; x<D_phi_S.w; x++){ 
	D_phi_S(x, y, t) = d_phi_S(sqr(G[0](x, y, 0)/hx) + sqr(G[0](x, y, 1)/hy) +
				   sqr(G[1](x, y, 0)/hx) + sqr(G[1](x, y, 1)/hy) +
				   sqr(G[2](x, y, 0)/hx) + sqr(G[2](x, y, 1)/hy));
      }		 
    }
  }
  return D_phi_S;
}


VariationalBasisProblemUV
VariationalBasisProblemUV::restriction(){
  int neww = (int)ceil(get_width()/2), newh = (int)ceil(get_height()/2);
  
  restrict = restriction_image(get_width(), get_height(), neww, newh);
    
  VariationalBasisProblemUV sub(neww, newh);
  sub.S_I.clear();
  sub.S_I_packers = S_I_packers;
  sub.times = times;
  sub.perFrameBeta = perFrameBeta;
  sub.basis = basis;

  for(int i=0; i<(int)S_I.size(); i++){
    sub.S_I.push_back(restrict.apply(S_I[i]));
  }
  sub.S_alpha = S_alpha;

  sub.D_0  = restrict.apply(D_0);
  sub.d    = restrict.apply(d);

  for (int i=0; i<(int)weights.size(); i++) {
    sub.weights.push_back(restrict.apply(weights[i]));
  }
  sub.hx =  hx * double(get_width())/double(neww);
  sub.hy =  hy * double(get_height())/double(newh);
  sub.beta = beta;
    
  sub.residuals = nacb::Imagef(neww, newh, residuals.nchannels);


  nacb::Imagef residuals_H;
  sub.f_h = nacb::Imagef(neww, newh, residuals.nchannels);
  sub.f_h = 0;
  sub.iteration(residuals_H, false);
    
  //residuals_H = 0 -AH*uv_H  => AH*uv_H = -residuals_H
  sub.f_h = restrict.apply(residuals) - residuals_H;// = r_H +  AH*uv_H;
  return sub;
}


robust_stats 
VariationalBasisProblemUV::iteration(nacb::Imagef & residuals, //Input output
				     bool update){
  nacb::StopWatch setup;

  nacb::Imagef d_old = d.copy();
  nacb::Imagef D = D_0 + d;
  std::vector<nacb::Imagef> D_grad(D.nchannels);
  std::vector<int> class_mapping;

  for (int i=0; i<D.nchannels; ++i) {
    D_grad[i] = D.gradient(i);
    class_mapping.push_back(i >= 0);
  }

  if (D.nchannels != 4) {
    static int warned = 0;
    if (warned == 0) {
      std::cerr << "Warning: wrong number of channels, fix mapping in D_phi_S!\n";
      std::cerr << "Warning: All basis functions (after the first will be robustified together.\n";
      warned = 1;
    }
  }
  int numUnknowns = D.nchannels;
  nacb::Imagef D_phi_S = get_D_phi_S(D_grad, hx, hy, d_phi_S, class_mapping, 2);

  nacb::Imagef D_phi_S_perFrame;
  if (perFrameBeta.size()) {
    D_phi_S_perFrame = get_D_phi_S_perFrame(D, hx, hy, basis, times, d_phi_S);
  }

  //printf("setup took %lf\n", double(setup));
  
  double res = 0.0;
  double len = 0.0;
  double moved = 0.0;

  int width = S_I[0].w, height = S_I[0].h;
  int nchans = S_alpha.size();
  
  if(residuals.w != width || residuals.h != height || residuals.nchannels != numUnknowns){
    residuals = nacb::Imagef(width, height, numUnknowns);
    residuals = 0;
  }

  nacb::Matrix basisMatrix = basis->getBasisMatrix(times);
#ifdef USE_OMP
#pragma omp parallel 
  {
#endif

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
  nacb::StopWatch timer;

#ifdef USE_OMP
#pragma omp parallel for
#endif


  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){	
      double d_phi_D_xy[nchans];
      
      // Current solution
      for (int i=0; i<D_0.nchannels; ++i)
	d_xy_vec[i] = d_old(x, y, i);
      d_xy_vec[D_0.nchannels] = 1;

      for(int si=0; si<nchans; si++){
	//unpackSymmetricMatrix<float>(&(S_I[si](x, y, 0)), S_I[si].nchannels, S[si]);
	S_I_packers[si].unpack(&(S_I[si](x, y, 0)), S_I[si].nchannels, S[si]);

	//S[si].printMatlab("SI");
	d_phi_D_xy[si] = d_phi_D(quadratic_form_eval(S[si], d_xy_vec, quad_form_temp));
      }

      int neigh[4][2] = {{-1,0}, {1,0}, {0,-1}, {0, 1}};
      double hs[4] = {hx, hx, hy, hy};
	
      for(int j=0; j<rhs.m; j++)
	rhs[j] = f_h(x, y, j);

      M.setAll(0);
	
      for(int si=0; si<(int)S_I.size(); si++){
	// !<per-pixel weighting function (on the data terms only).
	double wt_xy = weights[si](x, y);

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


      bool basisRegularization = !perFrameBeta.size();

      // Add a tiny bit of regularization for the case when we are performing per-frame flow
      // regularization.
      if (!basisRegularization)
	Msum.setAll(1e-6);

      for(int ni=0; ni<4; ni++){	
	int xn = x + neigh[ni][0];
	int yn = y + neigh[ni][1];
	
	// Reflect the boundaries.
	if(xn<0)xn = 0;
	if(yn<0)yn = 0;
	if(xn>=width)xn = width-1;
	if(yn>=height)yn = height-1;
	
	// Use most up to date solution	
	float * dd_xnyn = &d(xn, yn, 0);
	
	// Assume that the first element is a displacement (if the basis has one).
	// Otherwise, only regularize the channels if basisRegularization is true.
	int numBasisReg = basisRegularization ? D_0.nchannels : basis->hasDisplacementComponent();
	for(int j=0; j<numBasisReg; j++) {
	  int class_ind = j>=1;
	  double d_phi_S_xy = D_phi_S(x, y, class_ind);
	  double d_phi_S_xnyn = D_phi_S(xn, yn, class_ind);
	  
	  double d_phi_S_avg = (d_phi_S_xy + d_phi_S_xnyn)/2.0;
	  Msum[j] += beta.data[class_ind]*d_phi_S_avg/(hs[ni]*hs[ni]);
	  
	  // What would be used in the actual equations.
	  d_S_xy[j] += -beta.data[class_ind] * ((d_phi_S_xy + d_phi_S_xnyn)/2.0 * 
						(D(xn, yn, j) - D(x, y, j))/(hs[ni]*hs[ni]));
	  
	  
	  
	  rhs[j] += beta.data[class_ind]*d_phi_S_avg*(D_0(xn, yn, j) + dd_xnyn[j] - D_0(x, y, j))/(hs[ni]*hs[ni]);
	}
      }

      if (perFrameBeta.size()) {

	for(int ni=0; ni<4; ni++) {	
	  int xn = x + neigh[ni][0];
	  int yn = y + neigh[ni][1];
	  
	  // Reflect the boundaries.
	  if(xn<0)xn = 0;
	  if(yn<0)yn = 0;
	  if(xn>=width)xn = width-1;
	  if(yn>=height)yn = height-1;
	  
	  float * dd_xnyn = &d(xn, yn, 0);

	  // For each time step, we have a term -beta_ * phi(|\nabla u(t)|^2 + |\nabla v(t)|^2 + |\nabla w(t)|^2)
	  // In the euler lagrange, the term looks like -beta \nabla \dot (phi(*) b_i(t) (\nabla u(t)))
	  // Typically only one of u, v, t is affected by the basis.
	  for (int t=0; t<(int)times.size(); t++) {
	    double betaValue = perFrameBeta[t];
	    double d_phi_S_xy = D_phi_S_perFrame(x, y, t);
	    double d_phi_S_xnyn = D_phi_S_perFrame(xn, yn, t);		

	    int startChannel = 0; // Increase this if you don't want to do time-regularization of lower channels.
	    for (int i=startChannel; i<D_0.nchannels; i++) {
	      // Find which coordinate this channel influences (u, v, w)
	      // Get b_i(t)
	      // Basis matrix is numDims * numFrames*4
	      int influences = -1;
	      for (int k=1; k<4; k++) {
		if (fabs(basisMatrix(i, 4*t + k)) > 1e-5) {
		  influences = k;
		  break;
		}
	      }
	      if (influences < 0) 
		continue;
	      
	      double b_it = basisMatrix(i, 4*t + influences);
	      double d_phi_S_avg = (b_it*(d_phi_S_xy + d_phi_S_xnyn)/2.0);

	      for (int j=startChannel; j<D_0.nchannels; j++) {
		double b_jt = basisMatrix(j, 4*t + influences);
		// M(i, j) is non-zero iff it affects derivative at time t (e.g., if b_i(t) != 0)
		if (fabs(b_jt) < 1e-10) continue;

		M(i, j) += betaValue*d_phi_S_avg*b_jt/(hs[ni]*hs[ni]);
		rhs[i] += betaValue*d_phi_S_avg*b_jt*(D_0(xn, yn, j) + dd_xnyn[j] - D_0(x, y, j))/(hs[ni]*hs[ni]);

		// What would be used in the actual equations (used below for residual calculations)
		d_S_xy[i] += -betaValue * d_phi_S_avg * b_jt * (D(xn, yn, j) - D(x, y, j))/(hs[ni]*hs[ni]);
	      }
	    }
	  }
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
	double wt_xy = weights[si](x, y);
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

#ifdef USE_OMP
  } // End pragman omp
#endif

  // printf("%dx%d: %lf\n", f_h.w, f_h.h, double(timer));

  double erel = res/len;  
  robust_stats ret(moved, res, len, erel);

  //printf("%d %d\n", f_h.w, f_h.h);
  //ret.print();
  return ret;
}
