#include "math.h"
#include <boost/bind.hpp>
#include <boost/format.hpp>
//Yes, I am using two separate libraries for lbfgs, 
//this one is the packaged one (with bounds)
#include <lbfgs.h>

static uint8_t * bits16map = 0;

static void freeMap(){
  delete [] bits16map;
}

Mat4x4 getProjection4x4(Mat3x3 & K, Mat4x4 & E){
  Mat4x4 k4=Mat4x4::zeros();
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      k4(i,j)=K(i,j);
  Mat4x4 p=k4*E;
  for(int j=0;j<4;j++)
    p(3,j)=p(2,j);
  return p;
}

double xorFuncEvalAsym(const Image8 & psil, const Image8 & sils){
  double score = 0;
  if(bits16map==0){
    bits16map = new uint8_t[65536];
    
    for(int i=0; i<256; i++){
      for(int j=0; j<256; j++){
	int sij = 0;
	for(int k=0; k<8; k++){
	  int mask = 1<<k;
	  if( (i&mask) && !(j&mask)){
	    sij +=1;
	  }
	  else if ((j&mask) && !(i&mask)){
	    sij +=4;
	  }
	}
	bits16map[i*256+j] = sij;
      }
    }
  }
  int top=sils.w*sils.h*sils.nchannels; //Added the channels, shouldn't affect adversly
  for(int i=0;i<top;i++){
    int y = psil.data[i];
    int x = sils.data[i];
    int index = y*256+x;
    score+=bits16map[index];
  }
  return score;
}


double xorFuncEval8(const Image8 & psil, const Image8 & sils){
  static uint8_t bitcnt8[256];
  static int called=0;
  if(called==0){
    for(uint32_t c=0;c<256;c++){
      bitcnt8[c]=(c&0x1)+((c>>1)&0x1)+((c>>2)&0x1)+((c>>3)&0x1)+
	((c>>4)&0x1)+((c>>5)&0x1)+((c>>6)&0x1)+((c>>7)&0x1);
    }
    called=1;
  }
  int score=0;
  int top=sils.w*sils.h*sils.nchannels; //Added the channels, shouldn't affect adversly
  for(int i=0;i<top;i++){
    score+=bitcnt8[sils.data[i]^psil.data[i]];
    /*
      score+=(s)&0x1;
      score+=(s>>1)&0x1;
      score+=(s>>2)&0x1;
      score+=(s>>3)&0x1;
      score+=(s>>4)&0x1;
      score+=(s>>5)&0x1;
      score+=(s>>6)&0x1;
      score+=(s>>7)&0x1;
    */
  }
  /*
    for(int k=0;k<images.size();k++){
    for(int y=0;y<sils.h;y++){
    for(int x=0;x<sils.w;x++){
    score+=((sils(x,y)>>k)&0x1)!=(images[k](x,y,3)>128);
    sils(x,y)&psil(x,y);
    }
    }
    }
  */
  //printf("score is %d, t=%f (took %f)\n",score,t,(double)wtch);
  return score;
}


double jaccard(const Image8& psil, const Image8 & sil, int numViews) {
  double avg = 0;
  
  for(int i=0; i<numViews; i++){
    int is = 0;
    int un = 0;
    
    for(int y=0; y<psil.h; y++){
      for(int x=0; x<psil.w; x++){
	bool b1 = (psil(x, y)>>i)&0x1;
	bool b2 = (sil(x, y)>>i)&0x1;
	
	is += b1&&b2;
	un += b1||b2;
      }
    }
    avg += double(is)/double(un);
  }
  return double(avg)/numViews;
}


double lincomb(const boost::function<double (const Matrix &)> & fnc,
	       Matrix & p, Matrix & dir, double t){
  return fnc(p + dir*t);
}

Matrix powell(const boost::function<double (const Matrix &)> & fnc, Matrix & t, int maxits,
	      const OptimizerLogger & logger){
  int nparms = t.m*t.n;
  Matrix pback;
  Matrix parms;
  Matrix dir(nparms,1);
  
  pback = t.copy();
  parms = t.copy();

  vector<Matrix> extras;

  logger("Starting powell optimization.\n");

  for(int its=0; its<maxits; its++){
    Matrix orig=parms.copy();
    double tsum=0;

    for(int i=0; i<nparms; i++){
      dir.setAll(0);
      dir[i]=1;
      
      logger(OptimizerLogger::SubIterationLevel, 
	     (boost::format("Doing parameter %d.\n") % i).str());
      
      double t=golden_section(boost::bind(lincomb, fnc, parms, dir,_1),0.0,0.025);

      logger(OptimizerLogger::SubIterationLevel, 
	     (boost::format("done parameter %d (it: %d) %f\n") % i % its % t).str()),
      parms=parms+dir*t;
      tsum+=fabs(t);
    }
    for(int i=0; i<extras.size(); i++){
      double t=golden_section(boost::bind(lincomb, fnc, parms, extras[i], _1),0.0,0.025);
      printf("done xtra parameter %d (it:%d) %f\n",i,its,t);
      parms=parms+extras[i]*t;
    }
    if(tsum==0){
      logger(OptimizerLogger::OptimizationLevel, "Early break.\n");
      break;
    }
    dir=parms-orig;
    double t=golden_section(boost::bind(lincomb, fnc, parms, dir, _1),0.0,0.025);
    //boost::bind(&MultiCamera::xorFunc,this,psil,mesh,parms,dir,_1),0.0,0.025);
    parms=parms+dir*t;

    extras.push_back(dir.copy());

    logger(OptimizerLogger::IterationLevel, "Iteration complete.\n");
    if(!logger.iterationComplete())
      break;
  }
  return parms;
}

void powell(const boost::function<double (Armature *,Mesh *,Matrix &,Matrix & ,double)> & fnc,
	    Armature * arm,Mesh * mesh){
  int nparms = arm->getNParams();
  Matrix pback(nparms,1);
  Matrix parms;
  Matrix dir(nparms,1);
  arm->getParams(pback.data);
    
  parms=pback.copy();

  vector<Matrix> extras;

  for(int its=0; its<10; its++){
    Matrix orig=parms.copy();
    double tsum=0;
    for(int i=0;i<nparms;i++){
      dir.setAll(0);
      dir[i]=1;
      printf("doing parameter %d\n",i);
      double t=golden_section(boost::bind(fnc,arm,mesh,parms,dir,_1),0.0,0.025);
      printf("done parameter %d (it:%d) %f\n",i,its,t);
      parms=parms+dir*t;
      tsum+=fabs(t);
    }
    for(int i=0; i<extras.size(); i++){
      double t=golden_section(boost::bind(fnc,arm,mesh,parms,extras[i],_1),0.0,0.025);
      printf("done xtra parameter %d (it:%d) %f\n",i,its,t);
      parms=parms+extras[i]*t;
    }
    if(tsum==0){
      printf("early break !\n");
      break;
    }
    dir=parms-orig;
    double t=golden_section(boost::bind(fnc,arm,mesh,parms,dir,_1),0.0,0.025);//boost::bind(&MultiCamera::xorFunc,this,psil,mesh,parms,dir,_1),0.0,0.025);
    parms=parms+dir*t;

    extras.push_back(dir.copy());
  }
  arm->setParams(parms.data);
  arm->animate(*mesh);
}


void graddesc(const boost::function<double (Armature * ,Mesh *,Matrix &,Matrix & ,double)> & fnc,Armature * arm,Mesh * mesh){
  int nparms = arm->getNParams();
  Matrix pback(nparms,1);
  Matrix parms;
  Matrix dir(nparms,1);
  Matrix g(nparms,1);
  Matrix s(nparms,1);
  arm->getParams(pback.data);

  parms=pback.copy();

  dir.setAll(0);
  s.setAll(0);

  for(int its=0;its<30;its++){
    double cur=fnc(arm,mesh,parms,dir,0);
    printf("%d: %f\n",its,cur);

    double mx=0;
      
    g.setAll(0);

    for(int i=0;i<nparms;i++){
      dir[i]=1;
	
      double sp=0.0125/4;
      double up=fnc(arm,mesh,parms,dir,sp);
      while((up-cur)==0 && sp<0.1){
	sp*=2;
	up=fnc(arm,mesh,parms,dir,sp);
      }
      g[i]=(up-cur)/sp;
	
      dir[i]=0;
       
      mx=std::max(mx,fabs(g[i]));
    }
    g*=(1.0/mx);
    g*=-1;

    if(s.dot(s)<1e-8 || its%nparms==0)
      s=g.copy();
    else
      s=g-s*(s.dot(g)/(s.dot(s)));
    printf("%f\n",s.dot(g));
    //g.printMatlab("g");

    double alpha=4;
    double best=cur;
    double bestalpha=0;

    while(alpha>=1e-10){
      double val=fnc(arm,mesh,parms,s,alpha);
      //printf("a:%f %f\n",alpha,val);
      if(val<best){
	bestalpha=alpha;
	best=val;
	break;
      }
      alpha*=0.9;
    }
    if(bestalpha<=1e-10)break;
    parms=parms+s*bestalpha;
  }
  arm->setParams(parms.data);
  arm->animate(*mesh);
}


void lbfgs(const boost::function<double (Armature *,Mesh *,Matrix &,Matrix & ,double)> & fnc,Armature * arm,Mesh * mesh, double fdspc){
  int n=arm->getNParams();
  int nparms=n;
  int m=4;
  double F=1.0;
  Matrix g(n,1);
  Matrix W(n*(2*m+1)+2*m,1);
  int iprint[2]={1,3};
  int iflag=0;
  double eps=1e-7;
  double xtol=1e-7;
  int diagco=false;
  Matrix diag(n,1);

  Matrix pback(nparms,1);
  Matrix parms;
  Matrix dir(nparms,1);
    
  arm->getParams(pback.data);
  parms=pback.copy();
    
  dir.setAll(0);

  int its=0;

  while(iflag==0 || iflag==1){
    g.setAll(0);
    
    F=fnc(arm,mesh,parms,dir,0.0);
    for(int i=0;i<nparms;i++){
      dir[i]=1;
      
      double sp = fdspc; //Was 0.0125, changed to be param for other optimizations.
      double up = fnc(arm,mesh,parms,dir,sp);

      //Try some bigger spacing if the difference is zero.
      while((up-F)==0 && sp<fdspc*20){
	sp*=2;
	up=fnc(arm,mesh,parms,dir,sp);
      }
      g[i]=(up-F)/sp;
      
      dir[i]=0;
    }
    g.printMatlab("g");
      
    
    lbfgs_(n, m, parms.data, F, g.data, diagco, diag.data, iprint, eps, xtol, W.data, iflag);
    its++;

    if(its==200)break;

    if(iflag==0){
      printf("optimiziation complete!\n");
      break;
    }
    else if(iflag<0){
      printf("optimization error (%d)\n",iflag);
      break;
    }
  }

  arm->setParams(parms.data);
  arm->animate(*mesh);
}



/* \brief Helper evaluation function.
    
    A bit unfortunate that you can not pass matrices directly, and
    need to use the copying.
*/
static double lbfgs_boost_adapter(int n, double * x, double * g, const void * ptr){
  printf("Evalulating cost function.\n");
  const boost::function<double (Matrix &, Matrix &)> & func 
    = *((const boost::function<double (Matrix &, Matrix &)> *)ptr);
  Matrix xm(n, 1), grad(n, 1);
  memcpy(xm.data, x, n*sizeof(double));

  double ret = func(xm, grad);
  memcpy(g, grad.data, n*sizeof(double));
  return ret;
}


Matrix fd_gradient(const boost::function<double (Matrix & )> & func, Matrix & x, double dx){
  Matrix g(x.m, x.n);
  func_and_fd_gradient(func, x, dx, g);
}

double func_and_fd_gradient(const boost::function<double (Matrix & )> & func, Matrix & x, double dx, Matrix & g){
  if(g.m*g.n != x.m*x.n)g = Matrix(x.m, x.n);

  printf("Calling fd_gradient.\n");

  double cur = func(x);
    
  for(int i=0; i<x.m*x.n; i++){
    double back = x.data[i];
    x.data[i] += dx;
    g.data[i] = (func(x) - cur)/dx;
    x.data[i] = back;
  }
  return cur;
}

void lbfgs(const boost::function<double (Matrix &, Matrix & gradient)> & func, Matrix & x){
  printf("Allocting in lbfgs.\n");
  int * nbd = new int[x.m*x.n];
  memset(nbd, 0, sizeof(int)*x.m*x.n);
  Matrix lu = Matrix::zeros(x.m, x.n);

  printf("Calling lbfgs.\n");
  lbfgs(lbfgs_boost_adapter, x.m*x.n, x.data, lu.data, lu.data, nbd, &func);
  printf("Done calling lbfgs.\n");

  delete [] nbd;
}
