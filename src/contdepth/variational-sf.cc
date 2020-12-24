#include "variational-sf.h"
#include "Math.h"

using namespace std;
using namespace nacb;
using namespace mg;


template <bool UseDisparity> Vec2d VariationalSceneFlowProblem<UseDisparity>::defaultBeta(0.5, 0.5);
template <bool UseDisparity> double VariationalSceneFlowProblem<UseDisparity>::defaultAlpha = 1.0;


template <bool UseDisparity>
VariationalSceneFlowProblem<UseDisparity> VariationalSceneFlowProblem<UseDisparity>::restriction(){
  int neww = (int)ceil(get_width()/2), newh = (int)ceil(get_height()/2);
  
  restrict = restriction_image(get_width(), get_height(), neww, newh);
    
  VariationalSceneFlowProblem sub(neww, newh);
  sub.S_I.clear();

  for(int i=0; i<(int)S_I.size(); i++){
    sub.S_I.push_back(restrict.apply(S_I[i]));
  }
  sub.S_alpha = S_alpha;

  sub.D_0  = restrict.apply(D_0);
  sub.d    = restrict.apply(d);
  sub.weight = restrict.apply(weight);
  sub.hx =  hx * double(get_width())/double(neww);
  sub.hy =  hy * double(get_height())/double(newh);
  sub.beta = beta;
    
  sub.residuals = nacb::Imagef(neww, newh, 4);


  nacb::Imagef residuals_H;
  sub.f_h = nacb::Imagef(neww, newh, 4);
  sub.f_h = 0;
  sub.iteration(residuals_H, false);
    
  //residuals_H = 0 -AH*uv_H  => AH*uv_H = -residuals_H
  sub.f_h = restrict.apply(residuals) - residuals_H;// = r_H +  AH*uv_H;
  return sub;
}


template <bool UseDisparity>
robust_stats VariationalSceneFlowProblem<UseDisparity>::iteration(nacb::Imagef & residuals, //Input output
								  bool update){
  nacb::StopWatch setup;

  nacb::Imagef d_old = d.copy();
  nacb::Imagef D = D_0 + d;
  nacb::Imagef D_grad[4] = {D.gradient(0),
			    D.gradient(1),
			    D.gradient(2),
			    D.gradient(3)};
			    
  nacb::Imagef D_phi_S = get_D_phi_S(D_grad, hx, hy, d_phi_S);

  //printf("setup took %lf\n", double(setup));
  
  double res = 0.0;
  double len = 0.0;
  double moved = 0.0;

  int width = S_I[0].w, height = S_I[0].h;
  int nchans = S_alpha.size();

  if(residuals.w != width || residuals.h != height || residuals.nchannels != 4){
    residuals = nacb::Imagef(width, height, 4);
    residuals = 0;
  }

  nacb::StopWatch timer;

  for(int rb=0; rb<=1; ++rb){
#pragma omp parallel 
  {
  nacb::Matrix M(4, 4);
  nacb::Matrix rhs(4, 1);
  nacb::Matrix d_mat(4, 1);
  nacb::Matrix upd(4, 1);
  nacb::Matrix diff(4, 1);

#pragma omp for
  for(int y=0; y<height; y++){
    for(int x=rb; x<width; x+=2){	
      // Image data term
      Mat5x5 S[nchans];
      double d_phi_D_xy[nchans];
      Vec5d d_xy_vec(d_old(x, y, 0), d_old(x, y, 1), d_old(x, y, 2), d_old(x, y, 3), 1);
	
      for(int si=0; si<nchans; si++){
	S[si] = unpack_sym_mat5x5(S_I[si], x, y);
	d_phi_D_xy[si] = d_phi_D(d_xy_vec.dot(S[si]*d_xy_vec));
      }
	
      double d_S_xy[4] = {0.0, 0.0, 0.0, 0.0};

      int neigh[4][2] = {{-1,0}, {1,0}, {0,-1}, {0, 1}};
      double hs[4] = {hx, hx, hy, hy};
	
      for(int j=0; j<4; j++)
	rhs[j] = f_h(x, y, j);

      // !<per-pixel weighting function (on the data terms only).
      double wt_xy = weight(x, y);

      M.setAll(0);
	
      for(int si=0; si<(int)S_I.size(); si++){
	for(int i=0; i<4; i++){
	  for(int j=0; j<4; j++){
	    M(i, j) += S[si](i, j)*(d_phi_D_xy[si]*S_alpha[si]*wt_xy);
	  }
	  rhs[i] -=  S[si](i, 4)*(d_phi_D_xy[si]*S_alpha[si]*wt_xy);
	}
	// + alpha*d_phi_D2_xy*T(0, 0);	
      }

      double Msum[4] = {0.0, 0.0, 0.0, 0.0};

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
	for(int j=0; j<4; j++) {
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
      
      M(0, 0) += Msum[0];
      M(1, 1) += Msum[1];
      M(2, 2) += Msum[2];
      M(3, 3) += Msum[3];

      for(int j=0; j<4; j++)
	d_mat[j] = d(x, y, j);
      
      if(isNormal(M) && isNormal(rhs)) {
	// upd = Matrix::LlinLeastSq(M, rhs);
	solve_symmetric_4x4(M, rhs, upd);
      }
      else
	copy_matrix_no_alloc(upd, d_mat);

      if(!isNormal(M))
	copy_matrix_no_alloc(upd, d_mat);

      diff_matrix_no_alloc(diff, upd, d_mat);
            
      if(update){
	
	for(int j=0; j<4; j++){
	  if (fabs(upd[j]) > 1)
	    upd[j] = (upd[j] < 0)? -1: 1;
	  // Update solution, omega allows for SOR-type solution
	  d(x, y, j) = d(x, y, j)*(1.0-omega) + omega*upd[j];
	}
      }
      moved += diff.dot(diff);
      
      // These quantities would be used in the actual equations, just checking residual here.
      Vec5d dd_xy(d(x, y, 0), d(x, y, 1), d(x, y, 2), d(x, y, 3), 1);
      Vec5d d_T_xy;
      
      for(int si=0; si<nchans; si++)
	d_T_xy += (S[si]*dd_xy)*(S_alpha[si]*wt_xy*d_phi_D_xy[si]);
      
      Vec5d res_vec(f_h(x, y, 0) - (d_T_xy[0] + d_S_xy[0]),
		    f_h(x, y, 1) - (d_T_xy[1] + d_S_xy[1]),
		    f_h(x, y, 2) - (d_T_xy[2] + d_S_xy[2]),
		    f_h(x, y, 3) - (d_T_xy[3] + d_S_xy[3]), 0);
      
      // Should work also, except M has been modified by the solver...
      //double r = rhs - (M*d(x, y));
      // printf("%f %f\n", r, res_0);
      
      res += res_vec.dot(res_vec);
      len += upd.dot(upd);
      
      // These should probably use the updated coefficients, e.g., be recomputed after the loop,
      // when all the du, dv displacements have been updated.
      residuals(x, y, 0) = res_vec[0];
      residuals(x, y, 1) = res_vec[1];
      residuals(x, y, 2) = res_vec[2];
      residuals(x, y, 3) = res_vec[3];
    }
  }
  }
  } // red-black
  // printf("%dx%d: %lf\n", f_h.w, f_h.h, double(timer));

  double erel = res/len;  
  robust_stats ret(moved, res, len, erel);

  // printf("%d %d\n", f_h.w, f_h.h);
  // ret.print();
  return ret;
}



template class VariationalSceneFlowProblem<true>;
template class VariationalSceneFlowProblem<false>;
