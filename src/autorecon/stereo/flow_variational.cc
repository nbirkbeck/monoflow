/*
  I recently revisited the variational optic flow papers and
  realized that I may have misinterpreted the multi-resolution
  scheme for solving the flow equations.


  In the combined local and global approach (Lucas/Kanade Meets Horn/Schunk)
  there are some more details that say you warp the results from
  the coarser scale.

  Third revision: adding in the separate robustification
  and hopefully a fully approximate multi-grid scheme.

  Issues: better command line control, other robust issues (see class documentation).

*/
#include <stdio.h>
#include <nmath/vec3.h>
#include <nmath/matrix.h>
#include <nmath/sparsematrix.h>
#include <nimage/image.h>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <nmisc/timer.h>

#include "multigrid.h"
#include "flow.h"
#include "flow_variational.h"
#include <nmisc/commandline.h>

#include "../recon_globals.h"
#include "../recon_geometry.h"

using namespace std;
using namespace nacb;
using namespace mg;

double sqr(double d){
  return d*d;
}


/** \
 */
nacb::Imagef toGrayscaleIgnoreAlpha(const nacb::Imagef & im1){
  if(im1.nchannels>1)
    return (im1.getChannel(0) + im1.getChannel(1) + im1.getChannel(2))*0.33333;
  return im1;
}


/** \brief add an entry into the matrix if it is not already there (add to it if it is)
 */
void add(vector<std::unordered_map<int, double> > & mat,int i, int j, double val)
{
  assert(i<mat.size());
  
  if(mat[i].count(j))
    mat[i][j]+=val;
  else
    mat[i][j] = val;
}


/** \brief Aren't you sick of writing this yet?
    One dimensional gaussian filters (separable filtering).
 */
void gfilts(int hs, double std, Imagef & xf, Imagef & yf)
{
  double sum = 0;
  xf = Imagef(2*hs+1,1,1);
  yf = Imagef(1,2*hs+1,1);

  for(int i=0; i<=hs; i++){
    double ii = double(i)*i;
    double v = exp(-ii/(2.0*std*std));
    xf(hs+i,0) = v;
    xf(hs-i,0) = v;
    yf(0,hs+i) = v;
    yf(0,hs-i) = v;

    if(i==0) sum +=v;
    else sum +=2*v;      
  }
  xf *= (1.0/sum);
  yf *= (1.0/sum);
}

/** \brief the laplacian of the argument.
 */
Imagef lap(const Imagef & arg){
  Imagef res(arg.w, arg.h, 1);
  for(int y=0; y<arg.h; y++){
    for(int x=0; x<arg.w; x++){
      int xm1 = (x>0)?x-1:x+1;
      int ym1 = (y>0)?y-1:y+1;
      int xp1 = (x+1<arg.w)?x+1:x-1;
      int yp1 = (y+1<arg.h)?y+1:y-1;
      
      res(x,y) = arg(xm1,y)+arg(xp1,y)+arg(x,yp1)+arg(x,ym1)-4*arg(x,y);
    }
  }
  return res;
}

/** \brief Compute the variational optic flow (horn/schnuck) with local/global
    extension (no temporal smoothness).
  
    Estimate is the current estimate of the flow.  Assumes image im2 is
    already warped (this could change in the future).
    
    The only remaining issue is the solution takes a bit long using the
    umfsolve method.  Also, the smoothness, alpha, should be a parameter
    and may be dependent on the scale.  The best way to do the scale dependency
    is to use the pixel spacing when computing the laplacian and derivatices. 
    Assuming 1 at highest resolution should effectively solve the problem---if
    there is one.
*/
Imagef flow_variational_single(const Imagef & im1, const Imagef & im2,
			       const Imagef & estimate = Imagef()){
  using std::unordered_map;

  Imagef im1use = im1.boxFilter(1).boxFilter(1);
  Imagef im2use = im2.boxFilter(1).boxFilter(1);

  int m = 2*im1.w*im1.h;
  Matrix b(m, 1);
  double alpha = 0.003;
  
  Imagef g = im1.gradient(0);
  Imagef ft = (im2-im1);

  vector<unordered_map<int, double> > mat(m);

  Imagef xf,yf;
  
  gfilts(2, 1.0, xf, yf);

  printf("filtering\n");
  Imagef  ftt = (ft*ft).convolve(xf).convolve(yf);
  Imagef   ff = (g*g).convolve(xf).convolve(yf);
  Imagef  fxy = (g.getChannel(0)*g.getChannel(1)).convolve(xf).convolve(yf);
  Imagef  gxf = (g.getChannel(0)*ft).convolve(xf).convolve(yf);
  Imagef  gyf = (g.getChannel(1)*ft).convolve(xf).convolve(yf);
  printf("done filtering\n");

  //Want to minimize total laplacian, so take into account
  //the current estimate.
  Imagef lapest;
  if(estimate.w == im1.w && 
     estimate.h == im1.h)
    lapest = lap(estimate);
  else {
    lapest = Imagef(im1.w, im1.h, 1);
    lapest = 0;
  }
  for(int y=0; y<im1.h; y++){
    for(int x=0; x<im1.w; x++){
      int xp1 = ((x+1)<im1.w)?(x+1):(x-1);
      int xm1 = (x>0)?(x-1):(x+1);

      int yp1 = ((y+1)<im1.h)?(y+1):(y-1);
      int ym1 = (y>0)?(y-1):(y+1);
      

      int row = 2*(y*im1.w+x);
      int u = 2*(y*im1.w+x);
      int v = 2*(y*im1.w+x)+1;

      add(mat, row, 2*(yp1*im1.w+x), 1.0);
      add(mat, row, 2*(ym1*im1.w+x), 1.0);
      add(mat, row, 2*(y*im1.w+xp1), 1.0);
      add(mat, row, 2*(y*im1.w+xm1), 1.0);

      add(mat, row, u, -ff(x,y,0)/alpha - 4.0);//-g(x,y,0)*g(x,y,0)/alpha - 4.0);
      add(mat, row, v, -fxy(x,y)/alpha);//-g(x,y,0)*g(x,y,1)/alpha);
      b[row] = gxf(x,y)/alpha - lapest(x,y);//g(x,y,0)*ft(x,y)/alpha;
            
      add(mat, row+1, (2*(yp1*im1.w+x)+1), 1.0);
      add(mat, row+1, (2*(ym1*im1.w+x)+1), 1.0);
      add(mat, row+1, (2*(y*im1.w+xp1)+1), 1.0);
      add(mat, row+1, (2*(y*im1.w+xm1)+1), 1.0);
      
      add(mat, row+1, u, -fxy(x,y)/alpha);//-g(x,y,0)*g(x,y,1)/alpha);
      add(mat, row+1, v, -ff(x,y,1)/alpha - 4.0);//-g(x,y,1)*g(x,y,1)/alpha - 4.0);

      b[row+1] = gyf(x,y)/alpha - lapest(x,y);//g(x,y,1)*ft(x,y)/alpha;
    }
  }
  
  vector<vector<pair<int, double> > > cols;
  for(int i=0; i<m; i++){
    unordered_map<int,double>::iterator it = mat[i].begin();
    unordered_map<int,double>::iterator last = mat[i].end();
    vector<pair<int, double> > col;
    while(it!=last){
      col.push_back(pair<int,double>((*it).first, (*it).second));
      it++;
    }
    cols.push_back(col);
  }
  ft.write("/tmp/ft.png");
  printf("solving\n");
  SparseMatrix smat(m, cols);
  printf("after init matrix\n");
  Matrix res = SparseMatrix::umfsolve(smat,b);//smat.lsqr(b);
  printf("got solution\n");fflush(stdout);
  Imagef disp(im1.w, im1.h, 2);
  for(int y=0; y<im1.h; y++){
    for(int x=0; x<im1.w; x++){
      disp(x,y,0) = res[2*(y*im1.w+x)];
      disp(x,y,1) = res[2*(y*im1.w+x)+1];
    }
  }

  (im1.copy()).write("/tmp/warped0.png");
  (flow_warp_back(disp, im2)).write("/tmp/warped1.png");  
  return disp;
}


nacb::Mat3x3 outer_product(const Vec3d & x){
  Mat3x3 ret;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      ret(i, j) = x.data[i]*x.data[j];
    }
  }
  return ret;
}


double phi_1en6(double s2){
  const double eps2 = 1e-6;
  return sqrt(s2 + eps2);
}


double phi_1en2(double s2){
  const double eps2 = 1e-2;
  return sqrt(s2 + eps2);
}


double phi_1en4(double s2){
  const double eps2 = 1e-4;
  return sqrt(s2 + eps2);
}


double d_phi_const(double s2){
  return 1.0;
}


double d_phi_1en6(double s2){
  const double eps2 = 1e-6;
  return 0.5/sqrt(s2 + eps2);
}


double d_phi_1en4(double s2){
  const double eps2 = 1e-4;
  return 0.5/sqrt(s2 + eps2);
}


double d_phi_1en2(double s2){
  const double eps2 = 1e-2;
  return 0.5/sqrt(s2 + eps2);
}



Imagef get_D_phi_S(const nacb::Imagef & U_grad, const nacb::Imagef & V_grad, double hx, double hy,
		   double (*d_phi_S)(double)){
  nacb::Imagef D_phi_S(U_grad.w, U_grad.h, 1);

  for(int y=0; y<D_phi_S.h; y++){
    for(int x=0; x<D_phi_S.w; x++){ 
      D_phi_S(x, y) = 
	d_phi_S(sqr(U_grad(x, y, 0)/hx) + sqr(U_grad(x, y, 1)/hy) + 
		sqr(V_grad(x, y, 0)/hx) + sqr(V_grad(x, y, 1)/hy));
    }
  }
  return D_phi_S;
}


Imagef get_D_phi_S(const nacb::Imagef & D_grad, double hx, double hy, double (*d_phi_S)(double)){
  nacb::Imagef D_phi_S(D_grad.w, D_grad.h, 1);

  for(int y=0; y<D_phi_S.h; y++){
    for(int x=0; x<D_phi_S.w; x++){ 
      D_phi_S(x, y) = d_phi_S(sqr(D_grad(x, y, 0)/hx) + sqr(D_grad(x, y, 1)/hy));
    }
  }
  return D_phi_S;
}


Mat3x3 unpack_sym_mat3x3(const nacb::Imagef & M, int x, int y){
  Mat3x3 m;
  m(0, 0) = M(x, y, 0);
  m(1, 0) = m(0, 1) = M(x, y, 1);
  m(2, 0) = m(0, 2) = M(x, y, 2);
  
  m(1, 1) = M(x, y, 3);
  m(2, 1) = m(1, 2) = M(x, y, 4);
  m(2, 2) = M(x, y, 5);

  return m;
}


Imagef flow_variational(const Imagef & im1, const Imagef & im2){
  Imagef im1g, im2g;
  
  im1g = toGrayscaleIgnoreAlpha(im1);
  im2g = toGrayscaleIgnoreAlpha(im2);

  vector<Imagef> p1;
  vector<Imagef> p2;

  ///im1g = im1g.resize(im1g.w/2, im1g.h/2);
  ///im2g = im2g.resize(im2g.w/2, im2g.h/2);
  
  p1.push_back(im1g);
  p2.push_back(im2g);
  
  for(int i=0; i<4; i++){
    p1.push_back(p1[i].resize(p1[i].w/2,p1[i].h/2));
    p2.push_back(p2[i].resize(p2[i].w/2,p2[i].h/2));
  }
  Imagef est;
  for(int i=p1.size()-1; i>=0; i--){
    if(i==(int)p1.size()-1){
      est = flow_variational_single(p1[i], p2[i]);
    }
    else {
      printf("resizing %d %d  ==> %d %d\n", est.w, est.h, p1[i].w, p1[i].h);
      est = est.resize(p1[i].w, p1[i].h);
      est *= 2.0;
      p2[i] = flow_warp_back(est, p2[i]);
      Imagef est_upd = flow_variational_single(p1[i], p2[i], est);
      est += est_upd;
    }
  }
  return est;
}


/**
   Notes:  should be easy to make color channels (done). Just add an extra S_I (or compute sum like
           the gradient) and add the weight to S_alpha. x'*S1*x + x'*S2*x = x'*(S1 + S2)*x', so
	   you don't even need extra terms. (e.g., no more S_I)

	   -The color channels haven't been tested all that much.

   Issues: -when solving on image pyramid it doesn't seem to work when the downsampling
            uses the proper image spacing (instead it uses a 1).

	   -hopefully the per-pixel weights are interpreted correctly (and multigrid)
	   -fundamental multigrid seems to be working (still need to use a generic function).
	   -assumes the image2 is warped back.  This is not entirely correct; should use
	   the gradient through the warp (or compose the update properly).
	   
*/
class VariationalFlowProblem {
public:
  constexpr static double omega = 1.25;
  restriction_image restrict;

  double beta;  //!<Weight of smoothness term.
  double gamma; //!<Weight of f-matrix term.
  vector<double>   S_alpha; //!<Weight of data terms
  vector<nacb::Imagef> S_I; //!<Data term tensors

  nacb::Imagef UV_0, uv;
  nacb::Imagef f_h;
  nacb::Imagef residuals;

  nacb::Imagef Flines; //!<Fundamental matrix constraints.
  nacb::Imagef weight; //!<Per-pixel weighting.

  double              hx, hy;

  double (* d_phi_D)(double);
  double  (* d_phi_S)(double);


  VariationalFlowProblem(int width, int height){ 
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;    
  }

  VariationalFlowProblem(const nacb::Imagef & im1, 
			 const nacb::Imagef & im2,
			 const nacb::Imagef & _UV_0,
			 double _alpha = 1.0,
			 double _beta = 6.0/255.0,
			 double _hx = 1.0, double _hy = 1.0,
			 double _gamma = 0.0) : 
    beta(_beta), gamma(_gamma), UV_0(_UV_0), hx(_hx), hy(_hy) {
    
    d_phi_D = d_phi_1en6;
    d_phi_S  = d_phi_1en6;


    f_h = nacb::Imagef(im1.w, im1.h, 2);
    uv = nacb::Imagef(UV_0.w, UV_0.h, 2);
    residuals = nacb::Imagef(im1.w, im1.h, 2);
    weight = nacb::Imagef(im1.w, im1.h, 1);
    weight = 1;

    for(int i=0; i<2; i++){
      S_I.push_back(nacb::Imagef(im1.w, im1.h, 6));
      S_I.back() = 0; 

      S_alpha.push_back(1);
    }
    S_alpha[1] = _alpha;

    //A_h(x_h) = f_h, in our case, all the RHS is zero (to start with anyway)
    uv = 0;
    f_h = 0;

    int nchannels = std::min(im2.nchannels, 3);

    
    nacb::Imagef Ix[nchannels];
    nacb::Imagef Iy[nchannels];
    nacb::Imagef im1_grad[nchannels];

    for (int k=0; k<nchannels; k++) {
      nacb::Imagef im2_grad = im2.gradient(k);

      Ix[k] = im2_grad.getChannel(0)*(1.0/hx);
      Iy[k] = im2_grad.getChannel(1)*(1.0/hy);

      im1_grad[k] = im1.gradient(k);      
    }

    nacb::Imagef Iz = im2 - im1;
    nacb::Imagef Ixx[nchannels];
    nacb::Imagef Ixy[nchannels];
    nacb::Imagef Iyx[nchannels];
    nacb::Imagef Iyy[nchannels];
    nacb::Imagef Ixz[nchannels];
    nacb::Imagef Iyz[nchannels];
    
    for (int k=0; k<nchannels; k++) {
      nacb::Imagef Ix_grad = Ix[k].gradient(0);
      nacb::Imagef Iy_grad = Iy[k].gradient(0);

      Ixx[k] = Ix_grad.getChannel(0)*(1.0/hx);
      Ixy[k] = Ix_grad.getChannel(1)*(1.0/hy);
    
      Iyx[k] = Iy_grad.getChannel(0)*(1.0/hx);
      Iyy[k] = Iy_grad.getChannel(1)*(1.0/hy);
    
      Ixz[k] = Ix[k] - im1_grad[k].getChannel(0)*(1.0/hx);
      Iyz[k] = Iy[k] - im1_grad[k].getChannel(1)*(1.0/hy);
    
      
      //Make sure it is symmetric
      Ixy[k] = (Ixy[k] + Iyx[k])*0.5;
      Iyx[k] = Ixy[k];
    }

    int width = im1.w, height = im1.h;
    
    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	for (int k=0; k<nchannels; k++) {
	  Vec3d  Igrad(Ix[k](x, y), Iy[k](x, y), Iz(x, y, k));
	  Vec3d Igradx(Ixx[k](x, y), Ixy[k](x, y), Ixz[k](x, y));
	  Vec3d Igrady(Iyx[k](x, y), Iyy[k](x, y), Iyz[k](x, y));
	  
	  //Image data term
	  Mat3x3 S = outer_product(Igrad);
	  Mat3x3 T = outer_product(Igradx) + outer_product(Igrady);
	  
	  S_I[0](x, y, 0) += S(0, 0);
	  S_I[0](x, y, 1) += S(0, 1);
	  S_I[0](x, y, 2) += S(0, 2);
	  S_I[0](x, y, 3) += S(1, 1);
	  S_I[0](x, y, 4) += S(1, 2);
	  S_I[0](x, y, 5) += S(2, 2);
	  
	  S_I[1](x, y, 0) += T(0, 0);
	  S_I[1](x, y, 1) += T(0, 1);
	  S_I[1](x, y, 2) += T(0, 2);
	  S_I[1](x, y, 3) += T(1, 1);
	  S_I[1](x, y, 4) += T(1, 2);
	  S_I[1](x, y, 5) += T(2, 2);
	}
      }
    }

    Flines = Imagef(width, height, 3);
    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	Flines(x, y, 0) = sqrt(2);
	Flines(x, y, 1) = sqrt(2);
	Flines(x, y, 2) = 0;
      }
    }
    gamma = 0.0;
  }

  void setDataWeight(const nacb::Imagef & weights){
    weight = weights.copy();
  }

  void setFundamentalMatrix(const Matrix & F, double _gamma){
    gamma = _gamma;
    
    printf("Setting f-matrix.\n");
    Flines = Imagef(uv.w, uv.h, 3);

    for(int y=0; y<uv.h; y++){
      for(int x=0; x<uv.w; x++){
	Vec3d xp = F*Vec3d(x, y, 1);
	
	double len = xp.x*xp.x + xp.y*xp.y;
	xp *= 1.0/sqrt(len);
	/* //Transform the line so it can operate on displacements only.
	  a*(x + u) + b*(y + v) + c = 0;
	  a*x + a*u + b*y + b*v + c = 0;
	  a*x + b*y + (a*u + b*v + c) = 0
	  (a*u + b*v + (c +  a*x + b*y)) = 0
	*/
	xp.z += xp.dot(Vec3d(x, y, 0));
	
	Flines(x, y, 0) = xp.x;
	Flines(x, y, 1) = xp.y;
	Flines(x, y, 2) = xp.z;

	
      }
    }
    printf("done setting f-matrix.\n");
  }

  
  VariationalFlowProblem restriction(){
    int neww = (int)ceil(get_width()/2), newh = (int)ceil(get_height()/2);

    restrict = restriction_image(get_width(), get_height(), neww, newh);
    
    VariationalFlowProblem sub(neww, newh);
    sub.S_I.clear();

    for(int i=0; i<(int)S_I.size(); i++){
      sub.S_I.push_back(restrict.apply(S_I[i]));
    }
    sub.S_alpha = S_alpha;

    sub.UV_0  = restrict.apply(UV_0);
    sub.uv    = restrict.apply(uv);
    sub.weight = restrict.apply(weight);
    sub.hx =  hx * double(get_width())/double(neww);
    sub.hy =  hy * double(get_height())/double(newh);
    sub.beta = beta;
    sub.gamma = gamma;
    
    //Need to fix up these guys
    sub.Flines = restrict.apply(Flines);

    sub.residuals = nacb::Imagef(neww, newh, 2);

    

    nacb::Imagef residuals_H;    
    sub.f_h = nacb::Imagef(neww, newh, 2);
    sub.f_h = 0;
    sub.iteration(residuals_H, false);
    
    //residuals_H = 0 -AH*uv_H  => AH*uv_H = -residuals_H
    sub.f_h = restrict.apply(residuals) - residuals_H;// = r_H +  AH*uv_H;
    return sub;
  }

  nacb::Imagef prolong(const nacb::Imagef & img){
    return restrict.prolong(img);
  }

  void smooth(){
    iteration(residuals, true);
  }

  robust_stats iteration(nacb::Imagef & residuals, //Input output
			 bool update = true){
    Matrix M(2, 2);
    Matrix rhs(2, 1);
    
    nacb::StopWatch setup;

    nacb::Imagef uv_old = uv.copy();
    nacb::Imagef UV = UV_0 + uv;
    nacb::Imagef U_grad = UV.gradient(0);
    nacb::Imagef V_grad = UV.gradient(1);
    nacb::Imagef D_phi_S = get_D_phi_S(U_grad, V_grad, hx, hy, d_phi_S);

    //printf("setup took %lf\n", double(setup));
  
    double res = 0.0;
    double len = 0.0;
    double moved = 0.0;

    int width = S_I[0].w, height = S_I[0].h;
    int nchans = S_alpha.size();

    if(residuals.w != width || residuals.h != height || residuals.nchannels != 2){
      residuals = nacb::Imagef(width, height, 2);
      residuals = 0;
    }

    

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){	
	//Image data term
	Mat3x3 S[nchans];
	double d_phi_D_xy[nchans];
	Vec3d du_xy_vec(uv_old(x, y, 0), uv_old(x, y, 1), 1);

	for(int si=0; si<nchans; si++){
	  S[si] = unpack_sym_mat3x3(S_I[si], x, y);
	  d_phi_D_xy[si] = d_phi_D(du_xy_vec.dot(S[si]*du_xy_vec));
	}
	
	double d_S_u_xy = 0.0, d_S_v_xy = 0.0;
	double d_phi_S_xy = D_phi_S(x, y);
	
	int neigh[4][2] = {{-1,0}, {1,0}, {0,-1}, {0, 1}};
	double hs[4] = {hx, hx, hy, hy};

	M.setAll(0);
	rhs[0] = f_h(x, y, 0);
	rhs[1] = f_h(x, y, 1);

	//!<per-pixel weighting function (on the data terms only).
	double wt_xy = weight(x, y);
	
	for(int si=0; si<(int)S_I.size(); si++){
	  M(0, 0) += wt_xy*d_phi_D_xy[si]*S[si](0, 0)*S_alpha[si];// + alpha*d_phi_D2_xy*T(0, 0);
	  M(0, 1) += wt_xy*d_phi_D_xy[si]*S[si](0, 1)*S_alpha[si];// + alpha*d_phi_D2_xy*T(0, 1);
	  M(1, 1) += wt_xy*d_phi_D_xy[si]*S[si](1, 1)*S_alpha[si];// + alpha*d_phi_D2_xy*T(1, 1);
	  
	  rhs[0] -=  wt_xy*d_phi_D_xy[si]*S[si](0, 2)*S_alpha[si];
	  rhs[1] -=  wt_xy*d_phi_D_xy[si]*S[si](1, 2)*S_alpha[si];
	}
	M(1, 0) = M(0, 1);// S and T symmetric: d_phi_D1_xy*S(1, 0) + alpha*d_phi_D2_xy*T(1, 0);

	double Msum = 0.0;

	for(int ni=0; ni<4; ni++){	
	  int xn = x + neigh[ni][0];
	  int yn = y + neigh[ni][1];
	  
	  //Reflect the boundaries.
	  if(xn<0)xn = 1;
	  if(yn<0)yn = 1;
	  if(xn>=width)xn = width-2;
	  if(yn>=height)yn = height-2;
	  
	  double du_xnyn = uv(xn, yn, 0), dv_xnyn = uv(xn, yn, 1);
	  double d_phi_S_xnyn = D_phi_S(xn, yn);
	
	  //What would be used in the actual equations.
	  d_S_u_xy += -beta * ((d_phi_S_xy + d_phi_S_xnyn)/2.0 * 
			       (UV(xn, yn, 0) - UV(x, y, 0))/(hs[ni]*hs[ni]));
	  
	  d_S_v_xy += -beta * ((d_phi_S_xy + d_phi_S_xnyn)/2.0 * 
			       (UV(xn, yn, 1) - UV(x, y, 1))/(hs[ni]*hs[ni]));
	  
	  double d_phi_S_avg = (d_phi_S_xy + d_phi_S_xnyn)/2.0;
	  Msum += beta*d_phi_S_avg/(hs[ni]*hs[ni]);
	  
	  //Use most up to date solution.
	  rhs[0] += beta*d_phi_S_avg*(UV_0(xn, yn, 0) + du_xnyn - UV_0(x, y, 0))/(hs[ni]*hs[ni]);
	  rhs[1] += beta*d_phi_S_avg*(UV_0(xn, yn, 1) + dv_xnyn - UV_0(x, y, 1))/(hs[ni]*hs[ni]);
	}
      
	M(0, 0) += Msum;
	M(1, 1) += Msum;

	double d_F_u_xy = 0;
	double d_F_v_xy = 0;
	
	if(gamma>0){
	  Vec3d ln(Flines(x, y, 0), Flines(x, y, 1), Flines(x, y, 2));
	  double dot = ln.dot(Vec3d(UV(x, y, 0), UV(x, y, 1), 1.0));
	  double d_phi_F = d_phi_const(dot*dot);

	  M(0, 0) += wt_xy*gamma*d_phi_F*ln.x*ln.x;
	  M(0, 1) += wt_xy*gamma*d_phi_F*ln.x*ln.y;
	  M(1, 0) += wt_xy*gamma*d_phi_F*ln.y*ln.x;
	  M(1, 1) += wt_xy*gamma*d_phi_F*ln.y*ln.y;

	  rhs[0] -= wt_xy*gamma*d_phi_F*(ln.x*ln.x*UV_0(x,y,0) + ln.x*ln.y*UV_0(x,y,1) + ln.x*ln.z);
	  rhs[1] -= wt_xy*gamma*d_phi_F*(ln.y*ln.x*UV_0(x,y,0) + ln.y*ln.y*UV_0(x,y,1) + ln.y*ln.z);
	  
	  d_F_u_xy = wt_xy*gamma*d_phi_F*dot*ln.x;
	  d_F_v_xy = wt_xy*gamma*d_phi_F*dot*ln.y;
	}
	double det = (M(0,0)*M(1,1) - M(0,1)*M(0,1));

	//Matrix upd = Matrix::LlinLeastSq(M, rhs);
      
	double upd[2] = {( M(1,1)*rhs[0] - M(0,1)*rhs[1])/det,
			 (-M(0,1)*rhs[0] + M(0,0)*rhs[1])/det};
      
	Vec2d diff(upd[0] - uv(x, y, 0),
		   upd[1] - uv(x, y, 1));

	if(update){
	  //Update solution, omega allows for SOR-type solution
	  uv(x, y, 0) = uv(x, y, 0)*(1.0-omega) + omega*upd[0];
	  uv(x, y, 1) = uv(x, y, 1)*(1.0-omega) + omega*upd[1];
	}
	moved += diff.dot(diff);

	//These quantities would be used in the actual equations, just checking residual here.
	double du_xy = uv(x, y, 0), dv_xy = uv(x, y, 1);
	double d_T_u_xy = 0;
	double d_T_v_xy = 0;
	
	for(int si=0; si<nchans; si++){
	  d_T_u_xy += wt_xy*d_phi_D_xy[si]*(S[si](0,0)*du_xy + S[si](0,1)*dv_xy + S[si](0,2))*S_alpha[si];
	  d_T_v_xy += wt_xy*d_phi_D_xy[si]*(S[si](1,0)*du_xy + S[si](1,1)*dv_xy + S[si](1,2))*S_alpha[si];
	}
	double res_0 = f_h(x, y, 0) - (d_T_u_xy + d_S_u_xy + d_F_u_xy);
	double res_1 = f_h(x, y, 1) - (d_T_v_xy + d_S_v_xy + d_F_v_xy);
	
	//Should work also
	//double res_0 = rhs[0] - (M(0, 0)*uv(x, y, 0) + M(0, 1)*uv(x, y, 1));
	//double res_1 = rhs[1] - (M(1, 0)*uv(x, y, 0) + M(1, 1)*uv(x, y, 1));
	
	res += sqr(res_0) + sqr(res_1);
	len += sqr(Vec2d(upd[0], upd[1]).len());

	//These should probably use the updated coefficients, e.g., be recomputed after the loop,
	// when all the du, dv displacements have been updated.
	residuals(x, y, 0) = res_0;
	residuals(x, y, 1) = res_1;
      }
    }
    double erel = res/len;

    return robust_stats(moved, res, len, erel);
  }
  
  void add_to_solution(const nacb::Imagef & update){uv += update;}
  nacb::Imagef & get_solution(){return uv;}
  nacb::Imagef copy_solution(){return uv.copy();}

  int get_width() const {return uv.getWidth();}
  int get_height() const {return uv.getHeight();}
  bool is_linear() const {return false;}
  bool can_restrict() const {return get_width()>=16 && get_height()>=16;}


  nacb::Imagef getCombinedResult(){
    // Ideally, should probably be composed:
    // return flow_compose(UV_0, uv);
    // But it didn't help the results all that much, I think
    // it should still compute the gradients through the warped image.
    return UV_0 + uv;
  }
};



nacb::Imagef  flow_variational_robust(const nacb::Imagef & im1in,
				      const nacb::Imagef & im2in,
				      double alpha = 1.0,
				      double beta = 6.0/255.0, 
				      bool alphaMask = true,
				      int iters = 1){
  bool writing = true;
  mg::wcycle_solver solver;
  int outer_its = 1;
  solver.presmooth_its = solver.postsmooth_its = 5;
  solver.ngrids = 4;
  solver.wcycles = 1;

  nacb::Imagef weights;
  bool haveWeights = false;
  if(im1in.nchannels==4 && alphaMask) {
    weights = im1in.getChannel(3);
    haveWeights = true;
  }

  const bool grayScale = false;
  nacb::Imagef im1, im2;
  if (grayScale) {
    im1 = toGrayscaleIgnoreAlpha(im1in);
    im2 = toGrayscaleIgnoreAlpha(im2in);
  }
  else {
    im1 = im1in;
    im2 = im2in;
  }
    
  int ds = 16;
  nacb::Imagef disp(1, 1, 2);
  nacb::StopWatch totalTime;

  while(ds>=1){
    int w = im1.w/ds, h = im2.h/ds;
      
    nacb::Imagef im1_level = im1.resize(w, h);
    nacb::Imagef im2_level = im2.resize(w, h);
    nacb::Imagef wts_level;

    if(haveWeights)
      wts_level = weights.resize(w, h);

    if(disp.w==1){
      disp = nacb::Imagef(w, h, 2);
      disp = 0;
    }
    else {
      disp = disp.resize(w, h)*2.0;
    }

    for(int its=0; its<iters; its++){
      printf("on level %dx%d\n", w, h);
      nacb::Imagef im2w = flow_warp_back(disp, im2_level);
	
      VariationalFlowProblem prob(im1_level, im2w, disp, alpha, beta);//, ds, ds);

      if(haveWeights)prob.setDataWeight(wts_level);

      solver.solve(prob, outer_its);
      disp = prob.getCombinedResult();
    }

    if(writing){
      (im1_level).write("/tmp/i0.tga");
      (im2_level).write("/tmp/i1.tga");
      (flow_warp_back(disp, im2_level)).write("/tmp/w0.tga");
    }
    ds/=2;
  }
  printf("total time %lf\n", double(totalTime));
  return disp;
}

Matrix fundamental_euclidean(Matrix & F, const Matrix & K, 
			     Matrix & p1, Matrix & p2){
  Matrix U,S,V;

  Matrix W = Matrix::zeros(3, 3);
  Matrix Z = Matrix::zeros(3, 3);

  Z(0, 1) = 1;
  Z(1, 0) =-1; 

  W(0, 1) =-1;
  W(1, 0) = 1;
  W(2, 2) = 1;

  Matrix E = (K.transpose())*F*K;
  E.Lsvd(U, S, V);

  Matrix Rs[2] = {U*W*V.transpose(), 
		  U*W.transpose()*V.transpose()};
  Matrix ts[2] = {U.getColumn(2),
		  U.getColumn(2)*-1};
  Matrix Ext = Matrix::eye(4, 4);
  
  int best = 0;

  vector<Matrix> xs;
  xs.push_back(Matrix::cat(0, p1, Matrix::ones(1, p1.n)));
  xs.push_back(Matrix::cat(0, p2, Matrix::ones(1, p2.n)));
  

  //Two possibilities
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      //Matrix E = skewSymmetric(ts[j]) * Rs[i];
      printf("Tryign %d %d\n", i, j);
      vector<Matrix> Ps;
      Matrix E = getExtr(Rs[i], ts[j]);

      Ps.push_back(K*Matrix::eye(3, 4));
      Ps.push_back(K*Matrix::eye(3, 4)*E);

      printf("%dx%d\n", xs[0].m, xs[0].n);
      
      Matrix X = triangulate(Ps, xs);
      
      Matrix x0 = Ps[0]*X;
      Matrix x1 = Ps[1]*X;
      
      int infront = 0;
      
      for(int k=0; k<X.n; k++)
	infront += (x0(2, i)>0) && (x1(2, i)>0);
      
      printf("%d, %d  infront: %d/%d\n", i, j, infront, X.n);
      if(infront>=best){
	Ext = E;
	best = infront;

	if(infront == X.n)
	  Ext.printMatlab("Ext");

      }
    }
  }
  return Ext;
}


//FIXME: warning, this recomputes the grayscale images and the weight mask.
//FIXME: also does not take into account
nacb::Imagef flow_with_fundamental(const nacb::Imagef & im1in, 
				   const nacb::Imagef & im2in, 
				   const nacb::Imagef & dispin,
				   bool alphaMask = true){
  assert(dispin.nchannels == 2);

  mg::wcycle_solver solver;
  int outer_its = 1;
  solver.presmooth_its = solver.postsmooth_its = 5;
  solver.ngrids = 4;
  solver.wcycles = 1;

  nacb::Imagef weights;
  bool haveWeights = false;
  if(im1in.nchannels==4 && alphaMask) {
    weights = im1in.getChannel(3);
    haveWeights = true;
  }
  
  nacb::Imagef im1 = toGrayscaleIgnoreAlpha(im1in);
  nacb::Imagef im2 = toGrayscaleIgnoreAlpha(im2in);
  
  nacb::Imagef disp = dispin.copy();
  nacb::Image8 mask;

  if(haveWeights)
    mask = weights;
  else
    mask = (flow_length(disp)>2.0);  

  for(int after = 0; after<2; after++){
    Matrix p1, p2;
    Matrix F = flow_fundamental_matrix(disp, mask, &p1, &p2);
    F.printMatlab("F");

    Matrix K = Matrix::eye(3, 3);
    K(0, 0) = K(1, 1) = 366.73;
    K(0, 2) = 191.25;
    K(1, 2) = 149.54;
    
    fundamental_euclidean(F, K, p1, p2);
    
    nacb::Imagef im2w = flow_warp_back(disp, im2);
    VariationalFlowProblem prob(im1, im2w, disp, 1.0, 6.0/255.0);//, ds, ds);

    nacb::Imagef wts;
    double fwt = (after + 1);
    wts = (mask>0)*255;
    prob.setFundamentalMatrix(F, fwt);
    prob.setDataWeight(wts);

    solver.solve(prob, outer_its);

    disp = prob.getCombinedResult();

    (flow_warp_back(disp, im2)).write("/tmp/w0-f.tga");
    flow_save_length(disp, "/tmp/d-len-f.png");
  }
  mask.write("/tmp/mask.ppm");

  return disp;
}

#ifdef FLOW_VARIATIONAL_MAIN


int main(int ac, char * av[]){
  double grad = 1.0;
  double smooth = 6.0/255.0;
  int alphaMask = 0;
  int bidirectional = 0;
  int fmatrix = 0;
  int iters = 0;

  CommandLine cline;
  cline.registerOption("smooth", "Smoothness", &smooth, 's');
  cline.registerOption("grad", "Gradient weight", &grad, 'g');
  cline.registerOption("mask", "Alpha mask", &alphaMask, 0);
  cline.registerOption("bidi", "Bidirectional flow", &bidirectional, 0);
  cline.registerOption("fmatrix", "Fundamental matrix", &fmatrix, 0);
  cline.registerOption("iters", "Iterations", &iters, 0);
  cline.parse(ac, av);

  if(optind+1>=ac){
    printf("need two input images as arguments\n");
    return 0;
  }
  printf("Using images %s, %s\n", av[optind], av[optind+1]);

  Imagef im1(av[optind]);
  Imagef im2(av[optind+1]);

  if(1)
  {
    Imagef disp12 = flow_variational_robust(im1, im2, grad, smooth, alphaMask, iters);
    Imagef disp21;

    flow_save_flo(disp12, "/tmp/d12.flo");

    if(fmatrix)flow_with_fundamental(im1, im2, disp12, alphaMask);

    flow_save_length(disp12, "/tmp/d12-len.png");
    
    if(bidirectional){
      Imagef disp21 = flow_variational_robust(im2, im1, grad, smooth, alphaMask, iters);
      flow_save_length(disp21, "/tmp/d21-len.png");
      flow_save_flo(disp21, "/tmp/d21.flo");

      flow_save_crossfade("/tmp/cf-%04d.tga", im1, im2, disp12, disp21);
    }    
    return 0;
  }

  Imagef disp = flow_variational(im1, im2);
  flow_save_length(disp, "/tmp/d-len.png");
  return 0;
}


#endif
