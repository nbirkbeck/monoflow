#include "variational.h"

using namespace std;
using namespace nacb;
using namespace mg;


template <bool UseDisparity>
VariationalDispProblem<UseDisparity> VariationalDispProblem<UseDisparity>::restriction(){
  int neww = (int)ceil(get_width()/2), newh = (int)ceil(get_height()/2);
  
  restrict = restriction_image(get_width(), get_height(), neww, newh);
    
  VariationalDispProblem sub(neww, newh);
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
    
  sub.residuals = nacb::Imagef(neww, newh, 1);


  nacb::Imagef residuals_H;
  sub.f_h = nacb::Imagef(neww, newh, 1);
  sub.f_h = 0;
  sub.iteration(residuals_H, false);
    
  //residuals_H = 0 -AH*uv_H  => AH*uv_H = -residuals_H
  sub.f_h = restrict.apply(residuals) - residuals_H;// = r_H +  AH*uv_H;
  return sub;
}


template <bool UseDisparity>
robust_stats VariationalDispProblem<UseDisparity>::iteration(nacb::Imagef & residuals, //Input output
							     bool update){
  nacb::StopWatch setup;

  nacb::Imagef d_old = d.copy();
  nacb::Imagef D = D_0 + d;
  nacb::Imagef D_grad = D.gradient(0);
  nacb::Imagef D_phi_S = get_D_phi_S(D_grad, hx, hy, d_phi_S);

  //printf("setup took %lf\n", double(setup));
  
  double res = 0.0;
  double len = 0.0;
  double moved = 0.0;

  int width = S_I[0].w, height = S_I[0].h;
  int nchans = S_alpha.size();

  if(residuals.w != width || residuals.h != height || residuals.nchannels != 1){
    residuals = nacb::Imagef(width, height, 1);
    residuals = 0;
  }

  for(int rb=0; rb<=1; ++rb){
#pragma omp parallel for
  for(int y=0; y<height; y++){
    for(int x=rb; x<width; x+=2){	
      // Image data term
      Vec2d S[nchans];
      double d_phi_D_xy[nchans];
      Vec2d d_xy_vec(d_old(x, y), 1);
	
      for(int si=0; si<nchans; si++){
	S[si] = Vec2d(S_I[si](x, y, 0), S_I[si](x, y, 1));//unpack_sym_mat2x2(S_I[si], x, y);

	Vec2d Ssi2 = Vec2d(S_I[si](x, y, 1), S_I[si](x, y, 2));
	
	Vec2d Ssi_dx(S[si].dot(d_xy_vec), Ssi2.dot(d_xy_vec));
	d_phi_D_xy[si] = d_phi_D(d_xy_vec.dot(Ssi_dx));
      }
	
      double d_S_d_xy = 0.0;
      double d_phi_S_xy = D_phi_S(x, y);
	
      int neigh[4][2] = {{-1,0}, {1,0}, {0,-1}, {0, 1}};
      double hs[4] = {hx, hx, hy, hy};
	
      double rhs = f_h(x, y);

      //!<per-pixel weighting function (on the data terms only).
      double wt_xy = weight(x, y);
      double M = 0;
	
      for(int si=0; si<(int)S_I.size(); si++){
	M += wt_xy*d_phi_D_xy[si]*S[si].x*S_alpha[si];// + alpha*d_phi_D2_xy*T(0, 0);
	rhs -=  wt_xy*d_phi_D_xy[si]*S[si].y*S_alpha[si];
      }

      double Msum = 0.0;

      for(int ni=0; ni<4; ni++){	
	int xn = x + neigh[ni][0];
	int yn = y + neigh[ni][1];
	  
	// Reflect the boundaries.
	if(xn<0)xn = 0;
	if(yn<0)yn = 0;
	if(xn>=width)xn = width-1;
	if(yn>=height)yn = height-1;
	  
	double dd_xnyn = d(xn, yn);
	double d_phi_S_xnyn = D_phi_S(xn, yn);
	
	// What would be used in the actual equations.
	d_S_d_xy += -beta * ((d_phi_S_xy + d_phi_S_xnyn)/2.0 * 
			     (D(xn, yn) - D(x, y))/(hs[ni]*hs[ni]));
	  	  
	double d_phi_S_avg = (d_phi_S_xy + d_phi_S_xnyn)/2.0;
	Msum += beta*d_phi_S_avg/(hs[ni]*hs[ni]);
	  
	// Use most up to date solution.
	rhs += beta*d_phi_S_avg*(D_0(xn, yn) + dd_xnyn - D_0(x, y))/(hs[ni]*hs[ni]);
      }
      
      M += Msum;

      //M * upd = rhs
      //Matrix upd = Matrix::LlinLeastSq(M, rhs);

      if(M < 0 && M > -1e-7)
	M = -1e-7;
      else if(M > 0 && M < 1e-7)
	M = 1e-7;

      if(isnan(M) || isinf(M))
	M = 0;

      if(isnan(rhs) || isinf(rhs))
	rhs = 0;
      
      double upd = rhs/M;
      double diff = upd - d(x, y);

      // No change.
      if(isnan(upd) || isinf(upd))
	upd  = d(x, y);

      if(isnan(diff) || isinf(diff))
	diff = 0;

      if(update){
	// Update solution, omega allows for SOR-type solution
	d(x, y) = d(x, y)*(1.0-omega) + omega*upd;
      }
      moved += diff*diff;

      // These quantities would be used in the actual equations, just checking residual here.
      double dd_xy = d(x, y);
      double d_T_xy = 0;
	
      for(int si=0; si<nchans; si++){
	d_T_xy += wt_xy*d_phi_D_xy[si]*(S[si].x*dd_xy + S[si].y)*S_alpha[si];
      }
      double res_0 = f_h(x, y, 0) - (d_T_xy + d_S_d_xy);
	
      // Should work also
      // double r = rhs - (M*d(x, y));

      // printf("%f %f\n", r, res_0);

      res += res_0 * res_0;
      len += upd*upd;

      // These should probably use the updated coefficients, e.g., be recomputed after the loop,
      // when all the du, dv displacements have been updated.
      residuals(x, y) = res_0;
    }
  }
  } // red-black
  double erel = res/len;  
  return robust_stats(moved, res, len, erel);
}


template class VariationalDispProblem<true>;
template class VariationalDispProblem<false>;
