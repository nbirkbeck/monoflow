#include "recon_globals.h"
#include <assert.h>
#include <algorithm>
using namespace std;
using namespace nacb;


Vec2d distort(const Matrix & k,double a,double b){
  double r2= a*a+b*b;
  double r4 = r2*r2;
  double r6 = r4*r2;
  double a2=a*a;
  double b2=b*b;
  return Vec2d(a*(1.0+k[0]*r2+k[1]*r4 +k[4]*r6) + 2.0*k[2]*a*b + k[3]*(r2 + 2.0*a2),
	       b*(1.0+k[0]*r2+k[1]*r4 +k[4]*r6) + k[2]*(r2 + 2*b2) + 2*k[3]*a*b);
}
Vec2d intrinsic(const Matrix & K,
		const Vec2d & p){
  return Vec2d(K(0,0)*p.x+K(0,1)*p.y+K(0,2),
	                  K(1,1)*p.y+K(1,2));
}

template<class T>
Image<T> rectify(Image<T> & img,Matrix & K,Matrix & kc){
  Matrix Kinv = K.inverse();
  Image<T> out(img.w, img.h, img.nchannels);
  for(int y=0; y<img.h; y++){
    for(int x=0; x<img.w; x++){
      double xd = Kinv(0,0)*x+Kinv(0,1)*y+Kinv(0,2);
      double yd = Kinv(1,0)*x+Kinv(1,1)*y+Kinv(1,2);

      Vec2d pt = intrinsic(K, distort(kc, xd, yd));
      for(int k=0; k<img.nchannels; k++){
	out(x,y,k) = img.bilinear(pt.x, pt.y, k);
      }
    }
  }
  return out;
}

template Image8 rectify(Image8 & img,Matrix & K,Matrix & kc);
template Imagef rectify(Imagef & img,Matrix & K,Matrix & kc);


Image8 rgb2gray(Image8 & im)
{
  Image8 useim=im;
  if(im.nchannels>1){
    useim=Image8(im.w,im.h,1);
    for(int y=0;y<useim.h;y++){
      for(int x=0;x<useim.w;x++){
	useim(x,y)=(((short)im(x,y,0))+
		    ((short)im(x,y,1))+
		    ((short)im(x,y,2)))/3;
	
      }
    }
  }
  return useim;
}

void gradient(const Image8 & gim,Imagef & dx,Imagef & dy)
{/*
  Imagef gimf;
  gimf=gim;
  Imagef fx(3,1,1);
  Imagef fy(1,3,1);
  fx[0]=-1;
  fx[1]= 0;
  fx[2]= 1;
  fy[0]=-1;
  fy[1]= 0;
  fy[2]= 1;
  dx=gimf.convolve(fx);
  dy=gimf.convolve(fy);
*/
  
  dx=Imagef(gim.getWidth(),gim.getHeight(),1);
  dy=Imagef(gim.getWidth(),gim.getHeight(),1);

  for(int y=1;y<(gim.h-1);y++){
    for(int x=1;x<(gim.w-1);x++){
      dx(x,y)=(((float)gim(x+1,y))-((float)gim(x-1,y)))*0.5/255.0f;
      dy(x,y)=(((float)gim(x,y+1))-((float)gim(x,y-1)))*0.5/255.0f;
    }
  }
}

void gauss_sep(double sigma,Imagef & gx,Imagef & gy)
{
  int hw=(int)ceil(sigma*3);
  gx=Imagef(2*hw+1,1,1);
  gy=Imagef(1,2*hw+1,1);
  gx(hw,0)=1;
  gy(0,hw)=1;
  double sum=1;
  for(int i=0;i<hw;i++){
    double x=(i+1);
    gx(hw+i+1,0)=exp(-x*x/(2*sigma*sigma));
    gx(hw-i-1,0)=exp(-x*x/(2*sigma*sigma));
    gy(0,hw+i+1)=exp(-x*x/(2*sigma*sigma));
    gy(0,hw-i-1)=exp(-x*x/(2*sigma*sigma));
    sum+=2*exp(-x*x/(2*sigma*sigma));
  }
  gx*=(1.0/sum);
  gy*=(1.0/sum);
  /*
  for(int x=0;x<gx.w;x++)
    printf("%f ",gx(x,0));
  printf("\n");
  for(int y=0;y<gy.h;y++)
  printf("%f ",gy(0,y));*/
}

vector<Vec2f > harris(Image8 & im,int maxpts,int subpix)
{
  Image8 useim=im;
  if(im.nchannels>1){
    useim=Image8(im.w,im.h,1);
    for(int y=0;y<useim.h;y++){
      for(int x=0;x<useim.w;x++){
	useim(x,y)=(((short)im(x,y,0))+
		    ((short)im(x,y,1))+
		    ((short)im(x,y,2)))/3;
	
      }
    }
  }
  Imagef gx,gy;
  gradient(useim,gx,gy);

  //gx.writePNG("/tmp/gx.png");
  //gy.writePNG("/tmp/gy.png");

  float std=2;
  int hw=(int)ceil(3*std);
  int wsize=2*hw+1;
  Imagef gaussx(wsize,1,1);
  Imagef gaussy(1,wsize,1);
  
  float sum=1;
  gaussx(hw,0)=1;
  gaussy(0,hw)=1;
  for(int i=0;i<hw;i++){
    float x=(i+1)*(i+1);
    float val=exp(-x/(2*std*std));
    gaussx(i+hw+1,0)=val;
    gaussx(hw-i-1,0)=val;
    
    gaussy(0,i+hw+1)=val;
    gaussy(0,hw-i-1)=val;
    sum+=2*val;
  }
  gaussx*=(1.0/sum);
  gaussy*=(1.0/sum);
  
  /*printf("g=[");
  for(int i=0;i<gaussx.w;i++){
    printf("%f ",gaussx(i,0,0));
  }
  printf("];\n");*/

  Imagef gxx=(gx*gx).convolve(gaussx).convolve(gaussy);
  Imagef gxy=(gx*gy).convolve(gaussx).convolve(gaussy);
  Imagef gyy=(gy*gy).convolve(gaussx).convolve(gaussy);
  
  Imagef corn(useim.w,useim.h,1);
  float eps=1e-8;
  for(int y=0;y<corn.h;y++){
    for(int x=0;x<corn.w;x++){
      double det=gxx(x,y)*gyy(x,y)-gxy(x,y)*gxy(x,y);
      double tr=gxx(x,y)+gyy(x,y);
	   
      corn(x,y)=det/(tr+eps);

      if(isnan(corn(x,y)))corn(x,y)=0;
    }
  }
  //non-maximal supression (radius of 1)
  vector<PointAndScore> corns;
  for(int y=1;y<(corn.h-1);y++){
    for(int x=1;x<(corn.w-1);x++){
      if(corn(x,y)<corn(x-1,y))continue;
      if(corn(x,y)<corn(x+1,y))continue;
      if(corn(x,y)<corn(x,y+1))continue;
      if(corn(x,y)<corn(x,y-1))continue;
      if(corn(x,y)<corn(x+1,y+1))continue;
      if(corn(x,y)<corn(x+1,y-1))continue;
      if(corn(x,y)<corn(x-1,y+1))continue;
      if(corn(x,y)<corn(x-1,y-1))continue;
      if(corn(x,y)<2.0/65535.0f)continue;//discard bad points

      PointAndScore pad;
      pad.point=Vec2f(x,y);
      pad.score=-corn(x,y);
      corns.push_back(pad);
    }
  }
  //im.writePNG("/tmp/im.png");
  std::sort(corns.begin(),corns.end());
  printf("corns.size()=%ld\n", corns.size());
  printf("best: %f, worst: %f (%ld)\n",
         corns[0].score, corns[corns.size()-1].score, corns.size());
  vector<Vec2f > result;
  Matrix Ainv;
  Matrix b(9,1);
  if(subpix){
    Matrix A(9,6);
    int i=0;
    for(int y=-1;y<=1;y++){
      for(int x=-1;x<=1;x++){
	A(i,0)=x*x;A(i,1)=x*y;A(i,2)=y*y;
	A(i,3)=x;  A(i,4)=y;  A(i,5)=1;
	i++;
      }
    }
    Ainv=(A.transpose()*A).inverse()*A.transpose();
    //A.printMatlab("A");
    //Ainv.printMatlab("Ainv");
  }
  maxpts=std::min(maxpts, (int)corns.size());
  result.clear();
  for(int i=0; i < (int)corns.size() && (int)result.size() < maxpts; i++){
    result.push_back(corns[i].point);
    if(subpix){
      int cnt=0;
      for(int y=-1;y<=1;y++)
	for(int x=-1;x<=1;x++)
	  b[cnt++]=corn((int)round(corns[i].point.x)+x,(int)round(corns[i].point.y)+y);
      Matrix coeff=Ainv*b;
      //2ax+by+d=0
      //bx+2cy+e=0
      //x=(-d-by)/(2a)
      //y=(-2ea+bd)/(4*ac-b*b)

      double y0=(-2*coeff[4]*coeff[0]+coeff[1]*coeff[3])
	/(4*coeff[0]*coeff[2]-coeff[1]*coeff[1]);
      double x0=(-coeff[3]-coeff[1]*y0)/(2*coeff[0]);
      if(x0<=-1.0||y0<=-1.0||x0>=1.0||y0>=1.0){
	//printf("bad:%f,%f\n",x0,y0);
	x0=std::min(0.9999,std::max(x0,-0.9999));
	y0=std::min(0.9999,std::max(y0,-0.9999));
	result.pop_back();
	continue;
      }
      result[result.size()-1]+=Vec2f(x0,y0);
    }
  }
  //corn.writePNG("/tmp/corn.png");
  return result;
}


//FIXME: move this
template <class T>
Vec3<T> backProject(const Matrix & KRinv, const Matrix & Kt, T x, T y, T z){
  x*=z;
  y*=z;
  Matrix b(3,1);
  b[0]=x-Kt[0];
  b[1]=y-Kt[1];
  b[2]=z-Kt[2];
  Matrix X=KRinv*b;
  return Vec3f(X[0],X[1],X[2]);
}


template <class T>
Vec2<T> project(const Matrix & P,const Vec3<T> & pt,T * zout){
  double x,y,z;
  x=P.data[0]*pt.x+P.data[1]*pt.y+P.data[2]*pt.z+P.data[3];
  y=P.data[4]*pt.x+P.data[5]*pt.y+P.data[6]*pt.z+P.data[7];
  z=P.data[8]*pt.x+P.data[9]*pt.y+P.data[10]*pt.z+P.data[11];
  if(zout)*zout=z;
  return Vec2<T>(x/z,y/z);
}


Matrix projectPoints(const Matrix & P,const Matrix & X)
{
  Matrix p;
  if(X.m==4)p=P*X;
  else if(X.m==3){
    Matrix pts4(4,X.n);
    for(int j=0;j<X.n;j++){
      pts4(0,j)=X(0,j);
      pts4(1,j)=X(1,j);
      pts4(2,j)=X(2,j);
      pts4(3,j)=1;
    }
    p=P*pts4;
  }
  else throw std::string("X is wrong dimension in projectPoints");

  Matrix p2(2,p.n);
  for(int j=0;j<p.n;j++){
    p2(0,j)=p(0,j)/p(2,j);
    p2(1,j)=p(1,j)/p(2,j);
  }
  return p2;
}

/* Get the extrinsic matrix [R|t]
 */
Matrix getExtr(const Matrix & R,const Matrix & t)
{
  Matrix extr=Matrix::eye(4,4);
  extr(0,3)=t[0];
  extr(1,3)=t[1];
  extr(2,3)=t[2];
  for(int ii=0;ii<3;ii++){
    for(int jj=0;jj<3;jj++)
      extr(ii,jj)=R(ii,jj);
  }
  return extr;
}
/** Yet another projection.  This one has two terms
 *  for radial distortion...they may be 0.
 * \sa projectPoints
 */
Matrix projectPoints(const Matrix & A,const Matrix & dist,const Matrix & extr,const Matrix & pts)
{
  Matrix proj(2,pts.n);
  Matrix pts4(4,pts.n);
  for(int j=0;j<pts.n;j++){
    pts4(0,j)=pts(0,j);
    pts4(1,j)=pts(1,j);
    pts4(2,j)=pts(2,j);
    pts4(3,j)=1;
  }
  Matrix ptemp=Matrix::eye(3,4)*extr*pts4;
  for(int j=0;j<pts.n;j++){
    ptemp(0,j)=ptemp(0,j)/ptemp(2,j);
    ptemp(1,j)=ptemp(1,j)/ptemp(2,j);
    ptemp(2,j)=1;

    double r2=ptemp(0,j)*ptemp(0,j)+ptemp(1,j)*ptemp(1,j);
    double r4=r2*r2;
    double distortion=(dist[0]*r2+dist[1]*r4);
    ptemp(0,j)=ptemp(0,j)+ptemp(0,j)*distortion;
    ptemp(1,j)=ptemp(1,j)+ptemp(1,j)*distortion;
  }
  ptemp=A*ptemp;
  for(int j=0;j<pts.n;j++){
    proj(0,j)=ptemp(0,j);
    proj(1,j)=ptemp(1,j);
  }
  return proj;
}


/**
 *  Assuming that the given 2xnpts matrix represents an ideal 
 * projection (one with no distortion), we can compute the 
 * ideal normalized coords (which are those before projection).
 */
Matrix normalizeProjectedPoints(const Matrix & Ainv,const Matrix & pts)
{
  Matrix norm(2,pts.n);
  for(int j=0;j<pts.n;j++){
    norm(0,j)=Ainv(0,0)*pts(0,j)+Ainv(0,1)*pts(1,j)+Ainv(0,2);
    norm(1,j)=Ainv(1,0)*pts(0,j)+Ainv(1,1)*pts(1,j)+Ainv(1,2);
  }
  return norm;
}



/* \brief cost function for normalizeProjectedPoints
 */
static int norm_proj_cost(void * p,
			  int m,
			  int n,
			  const double * x,
			  double * fvec,
			  double * fjac,
			  int ldfjac,
			  int iflag){
  double * data = (double *) p;

  double r2 = x[0]*x[0] + x[1]*x[1];
  double r4 = r2*r2;
  double dd = (1.0 + r2*data[0] + r4*data[1]);
  double xd = x[0]*dd;
  double yd = x[1]*dd;

  if(iflag==1){
    fvec[0] = xd - data[2];
    fvec[1] = yd - data[3];
  }
  else if(iflag==2){
    double ddx = 2.0*data[0]*x[0] + 4.0*data[1]*r2*x[0];
    double ddy = 2.0*data[0]*x[1] + 4.0*data[1]*r2*x[1];

    fjac[0 +      0] = dd + x[0]*ddx;
    fjac[0 + ldfjac] = ddy;

    fjac[1 +      0] = ddx;
    fjac[1 + ldfjac] = dd + x[1]*ddy;
  }
  return iflag;
}

/** 
    \brief Compute normalized projected points when there is distortion.
    \param Ainv the inverse of the internal parameters
    \param d the distortion coefficients (only uses two terms)
*/
Matrix normalizeProjectedPoints(const Matrix & Ainv, 
				const Matrix & d,
				const Matrix & pts){
  Matrix fvec(2,1);
  Matrix xnorm = normalizeProjectedPoints(Ainv, pts);
  Matrix x(2,1);
  double data[4] = {d[0], d[1], 0, 0};

  for(int j=0; j<xnorm.n; j++){
    data[2] = xnorm(0,j);
    data[3] = xnorm(1,j);
    
    x[0] = xnorm(0,j);
    x[1] = xnorm(1,j);
        
    Matrix::lsqnonlin(norm_proj_cost, data, x, fvec);

    xnorm(0,j) = x[0];
    xnorm(1,j) = x[1];
  }
  return xnorm;
}

void factorProjectionMatrix(const Matrix & P, Matrix & K, Matrix & E){
  Matrix Q2,R2;
  Matrix perm(3,3);
  perm.setAll(0);
  perm(0,2) = perm(1,1) = perm(2,0) = 1;
  (P.submatrix(0,0,3,3).transpose()*perm).Lqr(Q2,R2);

  K = perm*R2.transpose()*perm;
  Matrix Q = perm*Q2.transpose();
  
  Matrix sc = Matrix::eye(3,3);
  if(K(0,0)<0)sc(0,0) = -1;
  if(K(1,1)<0)sc(1,1) = -1;
  if(K(2,2)<0)sc(2,2) = -1;
  
  K = K*sc;
  Q = sc*Q;
  
  E = Matrix::eye(4,4);
  Matrix t = K.inverse()*P.getColumn(3);
  E.setSubmatrix(0,0,Q);
  E.setSubmatrix(0,3,t);
}

double cubeeval(double a,double b,double c,double d,
		double x)
{
  return a*x*x*x+b*x*x+c*x+d;
}

int cuberoots(double a,double b,double c,double d,
	      double & x1,double & x2,double & x3)
{
  double f=(3.0*c/a-b*b/(a*a))/3.0;
  double g=(2.0*b*b*b/(a*a*a)-(9.0*b*c/(a*a))+27*d/a)/27.0;
  double h=g*g/4+f*f*f/27.0;

  if(f==0 && g==0 && h==0){
    double x1=-cbrt(d/a);
    x2=x1;
    x3=x1;
    return 1;//there is actually three, all the same
  }
  else if(h<=0){
    double i=sqrt((g*g/4) - h);
    double j=cbrt(i);
    double k=acos(-g/(2.0*i));
    double L=j*-1;
    double M=cos(k/3);
    double N=sqrt(3)*sin(k/3.0);
    double P=-(b/(3*a));
    x1=2*j*cos(k/3)-(b/(3.0*a));
    x2=L*(M+N)+P;
    x3=L*(M-N)+P;
    return 3;
  }
  //one real root
  else if(h>0){
    double R=-g/2.0+sqrt(h);
    double S=cbrt(R);
    double T=-g/2.0-sqrt(h);
    double U=cbrt(T);
    x1=(S+U)-b/(3.0*a);
    return 1;
  }
  return 0;
}

void unique_set(int * cset,int max,int n)
{
  assert(max>=n);
  for(int k=0;k<n;k++){
    int val=(int)(((double)rand())/((double)RAND_MAX)*(max));
    int unique=1;
    for(int k2=0;k2<k;k2++){
      if(val==cset[k2]){
	unique=0;
	break;
      }
    }
    if(!unique){
      k--;
      continue;
    }
    cset[k]=val;
  }
}

vector<int> unique(vector<int> & input)
{
  vector<int> result;
  for(int i=0;i<(int)input.size();i++){
    int good=1;
    for(int j=i+1;j<(int)input.size();j++){
      if(input[i]==input[j]){
	good=0;
	break;
      }
    }
    if(good)result.push_back(input[i]);
  }
  return result;
}



template Vec2f project(const Matrix & P,const Vec3f & pt,float * zout);
template Vec2d project(const Matrix & P,const Vec3d & pt,double * zout);


template Vec3f backProject(const Matrix & KRinv, const Matrix & Kt, float x, float y, float z);
template Vec3d backProject(const Matrix & KRinv, const Matrix & Kt, double x, double y, double z);
