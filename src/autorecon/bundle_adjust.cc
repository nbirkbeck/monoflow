#include "bundle_adjust.h"
#include "recon_globals.h"
#include <stdio.h>

using namespace std;

using namespace nacb;

namespace BA{
  int ncams;
  int npts;
  vector<Matrix> xs;
  double g_err;
  vector<vector<int> > obs;
  void pack(double * x,vector<Matrix> & Ps,Matrix & X){
    for(int i=0; i < (int)Ps.size();i++){
      memcpy(x+12*i,Ps[i].data,12*sizeof(double));
    }
    for(int i=0;i<X.n;i++){
      x[12*Ps.size()+3*i  ]=X(0,i);
      x[12*Ps.size()+3*i+1]=X(1,i);
      x[12*Ps.size()+3*i+2]=X(2,i);
    }
  }
  void unpack(double * x,vector<Matrix> & Ps,Matrix & X){
    X=Matrix(4,BA::npts);
    Ps.clear();
    for(int i=0;i<ncams;i++){
      Matrix P(3,4);
      memcpy(P.data,x+12*i,sizeof(double)*12);
      Ps.push_back(P);
    }
    for(int i=0;i<npts;i++){
      X(0,i)=x[12*ncams+3*i  ];
      X(1,i)=x[12*ncams+3*i+1];
      X(2,i)=x[12*ncams+3*i+2];
      X(3,i)=1;
    }
  }
  void bundleAdjustVis(const int & m,const int & n,
		       double * x,double * fvec,
		       double * fjac,const int & ldfjac,
		       int & info)
  {
    //unpack cameras and points
    Matrix X(4,BA::npts);
    vector<Matrix> Ps;
    unpack(x,Ps,X);
    
    //just get the error
    if(info==1){
      int ind=0;
      for(int i=0;i<ncams;i++){
	Matrix pp=projectPoints(Ps[i],X);
	for(int j=0;j<(int)obs[i].size();j++){
	  fvec[ind++]=pp(0,obs[i][j])-xs[i](0,j);
	  fvec[ind++]=pp(1,obs[i][j])-xs[i](1,j);
	}
      }
      g_err=0;
      for(int i=0;i<ind;i++){
	double d=fvec[i]*fvec[i];
	if(isinf(d)||isnan(d)){
	  d=1000000;
	  fvec[i]=d;
	  info=-1;
	}
	g_err+=d;
      }
      g_err/=(m);
    }
    //just update the Jacobian
    else if(info==2){
      for(int j=0;j<n;j++)
	for(int i=0;i<m;i++)
	  fjac[i+ldfjac*j]=0;
      int row=0;
      for(int i=0;i<ncams;i++){
	int col=12*i;

	Matrix PsX=Ps[i]*X;
	for(int j=0;j<(int)obs[i].size();j++){
	  int jind=obs[i][j];
	  double imz=PsX(2,jind);
	  double imzz=imz*imz;

	  fjac[    col*ldfjac+(row+2*j)]=X(0,jind)/imz;
	  fjac[(col+1)*ldfjac+(row+2*j)]=X(1,jind)/imz;
	  fjac[(col+2)*ldfjac+(row+2*j)]=X(2,jind)/imz;
	  fjac[(col+3)*ldfjac+(row+2*j)]=X(3,jind)/imz;
	  fjac[(col  )*ldfjac+(row+2*j+1)]=0;
	  fjac[(col+1)*ldfjac+(row+2*j+1)]=0;
	  fjac[(col+2)*ldfjac+(row+2*j+1)]=0;
	  fjac[(col+3)*ldfjac+(row+2*j+1)]=0;

	  fjac[(col+4)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+5)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+6)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+7)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+4)*ldfjac+(row+2*j+1)]=X(0,jind)/imz;
	  fjac[(col+5)*ldfjac+(row+2*j+1)]=X(1,jind)/imz;
	  fjac[(col+6)*ldfjac+(row+2*j+1)]=X(2,jind)/imz;
	  fjac[(col+7)*ldfjac+(row+2*j+1)]=X(3,jind)/imz;

	  fjac[ (col+8)*ldfjac+(row+2*j)]  =-X(0,jind)*PsX(0,jind)/imzz;
	  fjac[ (col+8)*ldfjac+(row+2*j+1)]=-X(0,jind)*PsX(1,jind)/imzz;

	  fjac[ (col+9)*ldfjac+(row+2*j)]  =-X(1,jind)*PsX(0,jind)/imzz;
	  fjac[ (col+9)*ldfjac+(row+2*j+1)]=-X(1,jind)*PsX(1,jind)/imzz;

	  fjac[(col+10)*ldfjac+(row+2*j)]  =-X(2,jind)*PsX(0,jind)/imzz;
	  fjac[(col+10)*ldfjac+(row+2*j+1)]=-X(2,jind)*PsX(1,jind)/imzz;

	  fjac[(col+11)*ldfjac+(row+2*j)]  =-X(3,jind)*PsX(0,jind)/imzz;
	  fjac[(col+11)*ldfjac+(row+2*j+1)]=-X(3,jind)*PsX(1,jind)/imzz;

	  int ptcol=jind*3+ncams*12;
	  fjac[    ptcol*ldfjac+(row+2*j  )]=
	    Ps[i](0,0)/PsX(2,jind)-Ps[i](2,0)*PsX(0,jind)/imzz;
	  fjac[    ptcol*ldfjac+(row+2*j+1)]=
	    Ps[i](1,0)/PsX(2,jind)-Ps[i](2,0)*PsX(1,jind)/imzz;

	  fjac[(ptcol+1)*ldfjac+(row+2*j  )]=
	    Ps[i](0,1)/PsX(2,jind)-Ps[i](2,1)*PsX(0,jind)/imzz;
	  fjac[(ptcol+1)*ldfjac+(row+2*j+1)]=
	    Ps[i](1,1)/PsX(2,jind)-Ps[i](2,1)*PsX(1,jind)/imzz;

	  fjac[(ptcol+2)*ldfjac+(row+2*j  )]=
	    Ps[i](0,2)/PsX(2,jind)-Ps[i](2,2)*PsX(0,jind)/imzz;
	  fjac[(ptcol+2)*ldfjac+(row+2*j+1)]=
	    Ps[i](1,2)/PsX(2,jind)-Ps[i](2,2)*PsX(1,jind)/imzz;
	}
	row+=2*obs[i].size();
      }
    }
  }
  void bundleAdjustErrAn(const int & m,const int & n,
		       double * x,double * fvec,
		       double * fjac,const int & ldfjac,
		       int & info)
  {
    //unpack cameras and points
    Matrix X(4,BA::npts);
    vector<Matrix> Ps;
    unpack(x,Ps,X);
    
    //just get the error
    if(info==1){
      for(int i=0;i<ncams;i++){
	Matrix pp=projectPoints(Ps[i],X);
	for(int j=0;j<npts;j++){
	  fvec[i*npts*2+2*j  ]=pp(0,j)-xs[i](0,j);
	  fvec[i*npts*2+2*j+1]=pp(1,j)-xs[i](1,j);
	}
      }
      g_err=0;
      for(int i=0;i<ncams*npts*2;i++){
	double d=fvec[i]*fvec[i];
	if(isinf(d)||isnan(d)){
	  d=1000000;
	  fvec[i]=d;
	  info=-1;
	}
	g_err+=d;
      }
      g_err/=(m);
      static int cnt=0;
      //if(cnt%10==0)printf("cnt:%d,%f\n",cnt,g_err);
      cnt++;
    }
    //just update the Jacobian
    else if(info==2){
      for(int j=0;j<n;j++)
	for(int i=0;i<m;i++)
	  fjac[i+ldfjac*j]=0;
     
      for(int i=0;i<ncams;i++){
	int col=12*i;
	int row=i*npts*2;
	Matrix PsX=Ps[i]*X;
	for(int j=0;j<npts;j++){
	  double imz=PsX(2,j);
	  double imzz=imz*imz;

	  fjac[    col*ldfjac+(row+2*j)]=X(0,j)/imz;
	  fjac[(col+1)*ldfjac+(row+2*j)]=X(1,j)/imz;
	  fjac[(col+2)*ldfjac+(row+2*j)]=X(2,j)/imz;
	  fjac[(col+3)*ldfjac+(row+2*j)]=X(3,j)/imz;
	  fjac[(col  )*ldfjac+(row+2*j+1)]=0;
	  fjac[(col+1)*ldfjac+(row+2*j+1)]=0;
	  fjac[(col+2)*ldfjac+(row+2*j+1)]=0;
	  fjac[(col+3)*ldfjac+(row+2*j+1)]=0;

	  fjac[(col+4)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+5)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+6)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+7)*ldfjac+(row+2*j)  ]=0;
	  fjac[(col+4)*ldfjac+(row+2*j+1)]=X(0,j)/imz;
	  fjac[(col+5)*ldfjac+(row+2*j+1)]=X(1,j)/imz;
	  fjac[(col+6)*ldfjac+(row+2*j+1)]=X(2,j)/imz;
	  fjac[(col+7)*ldfjac+(row+2*j+1)]=X(3,j)/imz;

	  fjac[ (col+8)*ldfjac+(row+2*j)]  =-X(0,j)*PsX(0,j)/imzz;
	  fjac[ (col+8)*ldfjac+(row+2*j+1)]=-X(0,j)*PsX(1,j)/imzz;

	  fjac[ (col+9)*ldfjac+(row+2*j)]  =-X(1,j)*PsX(0,j)/imzz;
	  fjac[ (col+9)*ldfjac+(row+2*j+1)]=-X(1,j)*PsX(1,j)/imzz;

	  fjac[(col+10)*ldfjac+(row+2*j)]  =-X(2,j)*PsX(0,j)/imzz;
	  fjac[(col+10)*ldfjac+(row+2*j+1)]=-X(2,j)*PsX(1,j)/imzz;

	  fjac[(col+11)*ldfjac+(row+2*j)]  =-X(3,j)*PsX(0,j)/imzz;
	  fjac[(col+11)*ldfjac+(row+2*j+1)]=-X(3,j)*PsX(1,j)/imzz;

	  int ptcol=j*3+ncams*12;
	  fjac[    ptcol*ldfjac+(row+2*j  )]=
	    Ps[i](0,0)/PsX(2,j)-Ps[i](2,0)*PsX(0,j)/imzz;
	  fjac[    ptcol*ldfjac+(row+2*j+1)]=
	    Ps[i](1,0)/PsX(2,j)-Ps[i](2,0)*PsX(1,j)/imzz;

	  fjac[(ptcol+1)*ldfjac+(row+2*j  )]=
	    Ps[i](0,1)/PsX(2,j)-Ps[i](2,1)*PsX(0,j)/imzz;
	  fjac[(ptcol+1)*ldfjac+(row+2*j+1)]=
	    Ps[i](1,1)/PsX(2,j)-Ps[i](2,1)*PsX(1,j)/imzz;

	  fjac[(ptcol+2)*ldfjac+(row+2*j  )]=
	    Ps[i](0,2)/PsX(2,j)-Ps[i](2,2)*PsX(0,j)/imzz;
	  fjac[(ptcol+2)*ldfjac+(row+2*j+1)]=
	    Ps[i](1,2)/PsX(2,j)-Ps[i](2,2)*PsX(1,j)/imzz;
	}
      }
    }
  }
  void bundleAdjustErr(const int & m,const int & n,
		       double * x,double * fvec,
		       int & info)
  {
    //unpack cameras and points
    Matrix X(4,BA::npts);
    vector<Matrix> Ps;
    unpack(x,Ps,X);

    for(int i=0;i<ncams;i++){
      Matrix pp=projectPoints(Ps[i],X);
      for(int j=0;j<npts;j++){
	fvec[i*npts*2+2*j  ]=pp(0,j)-xs[i](0,j);
	fvec[i*npts*2+2*j+1]=pp(1,j)-xs[i](1,j);
      }
    }
    g_err=0;
    for(int i=0;i<ncams*npts*2;i++){
      double d=fvec[i]*fvec[i];
      if(isinf(d)||isnan(d)){
	info=-1;
	d=1000000;
	fvec[i]=d;
      }
      g_err+=d;
    }
    g_err/=(m);
  }
}

#include <sba.h>

class EuclideanBACache{
public:
  Matrix * Rs;
  Matrix * R_axs;
  Matrix * R_ays;
  Matrix * R_azs;
  vector<Vec3d> keys;
  vector<Vec3d> keys_d;
  int hits[2];
  int misses[2];
  int dofocal;//not really a cache, but are we adjusting focal lengths?
  vector<Matrix> Ks;
  EuclideanBACache(int ncams=0,const vector<Matrix> & Ks=vector<Matrix>(),
		   int dofocal=0){
    this->Ks=Ks;
    //make sure that we have enough intrinsics
    for(int i=0;i<ncams;i++){
      if(i >= (int)this->Ks.size())
	this->Ks.push_back(Ks[0]);
    }
    for(int i=0;i<ncams;i++){
      keys.push_back(Vec3d(100000,-100000,100000));
      keys_d.push_back(Vec3d(100000,-100000,100000));
    }
    R_axs=new Matrix[ncams];
    R_ays=new Matrix[ncams];
    R_azs=new Matrix[ncams];
    Rs=new Matrix[ncams];

    hits[0]=0;
    hits[1]=0;
    misses[0]=misses[1]=0;
    this->dofocal=dofocal;
  }
  ~EuclideanBACache(){
    //printf("hits/accesses: %d/%d  %d/%d  %d\n",hits[0],hits[0]+misses[0],hits[1],hits[1]+misses[1],keys.size());
    delete [] R_axs;
    delete [] R_ays;
    delete [] R_azs;
    delete [] Rs;
  }

  Matrix & getAndCacheR(int cam,double * aj){
    if(keys[cam].x==aj[0] &&
       keys[cam].y==aj[1] &&
       keys[cam].z==aj[2]){
      hits[0]++;
    }
    else{
      Matrix rod(3,1);
      rod[0]=aj[0];
      rod[1]=aj[1];
      rod[2]=aj[2];
      Rs[cam]=rod.rod2matrix();
      keys[cam]=Vec3d(aj[0],aj[1],aj[2]);
      misses[0]++;
    }
    return Rs[cam];
  }
  void cacheDeriv(int cam,double * aj){
    getAndCacheR(cam,aj);
    if(keys_d[cam].x==aj[0] &&
       keys_d[cam].y==aj[1] &&
       keys_d[cam].z==aj[2]){
      hits[1]++;
    }
    else{
      double angle=sqrt(aj[0]*aj[0]+aj[1]*aj[1]+aj[2]*aj[2]);
      Matrix S(3,3);
      S(0,0)=0;
      S(0,1)=-aj[2]/angle;
      S(0,2)= aj[1]/angle;
      S(1,0)= aj[2]/angle;
      S(1,1)= 0;
      S(1,2)=-aj[0]/angle;
      S(2,0)=-aj[1]/angle;
      S(2,1)= aj[0]/angle;
      S(2,2)=0;
      Matrix SS=S*S;
      
      double sin_angle=sin(angle);
      double cos_angle=cos(angle);
      double angle_ax=aj[0]/angle;
      double angle_ay=aj[1]/angle;
      double angle_az=aj[2]/angle;
      Matrix S_ax(3,3);
      Matrix S_ay(3,3);
      Matrix S_az(3,3);
      Matrix R_ax;
      Matrix R_ay;
      Matrix R_az;
      
      S_ax.setAll(0);
      S_ay.setAll(0);
      S_az.setAll(0);
      S_ax(1,2)=-1;
      S_ax(2,1)= 1;
      S_ay(0,2)=1;
      S_ay(2,0)=-1;
      S_az(0,1)=-1;
      S_az(1,0)=1;
      S_ax=S_ax*(1.0/angle)-S*(angle_ax/angle);
      S_ay=S_ay*(1.0/angle)-S*(angle_ay/angle);
      S_az=S_az*(1.0/angle)-S*(angle_az/angle);
      
      R_ax=S_ax*sin_angle+S*(cos_angle*angle_ax)+
	(S_ax*S+S*S_ax)*(1.0-cos_angle)+SS*sin_angle*angle_ax;
      R_ay=S_ay*sin_angle+S*(cos_angle*angle_ay)+
	(S_ay*S+S*S_ay)*(1.0-cos_angle)+SS*sin_angle*angle_ay;
      R_az=S_az*sin_angle+S*(cos_angle*angle_az)+
	(S_az*S+S*S_az)*(1.0-cos_angle)+SS*sin_angle*angle_az;

      R_axs[cam]=R_ax;
      R_ays[cam]=R_ay;
      R_azs[cam]=R_az;
      keys_d[cam]=Vec3d(aj[0],aj[1],aj[2]);
      misses[1]++;
    }
  }
  
};

void eucproj(int j, int i,
	     double *aj, double *bi,
	     double *xij,void *adata){
  EuclideanBACache * cache=(EuclideanBACache *) adata;
  Matrix & R=cache->getAndCacheR(j,aj);
  const Matrix & K=cache->Ks[j];
  const double & px=K.get(0,2);
  const double & py=K.get(1,2);
  const double & fx=K.get(0,0);
  const double & fy=K.get(1,1);
  
    /*Matrix rod(3,1);
  rod[0]=aj[0];
  rod[1]=aj[1];
  rod[2]=aj[2];
  Matrix R=rod.rod2matrix();*/
  const double * Rd=R.data;
  double scf=1.0;
  if(cache->dofocal)scf=aj[6];

  double x= Rd[0]*bi[0]+Rd[1]*bi[1]+Rd[2]*bi[2]+aj[3];
  double y= Rd[3]*bi[0]+Rd[4]*bi[1]+Rd[5]*bi[2]+aj[4];
  double z= Rd[6]*bi[0]+Rd[7]*bi[1]+Rd[8]*bi[2]+aj[5];
  xij[0]=scf*fx*x/z+px;
  xij[1]=scf*fy*y/z+py;
}
void eucprojac(int j, int i,
	    double *aj, double *bi,
  double *Aij, double *Bij, void *adata){
  
  EuclideanBACache * cache=(EuclideanBACache *) adata;
  cache->cacheDeriv(j,aj);
  Matrix  R=cache->Rs[j];
  Matrix  R_ax=cache->R_axs[j];
  Matrix  R_ay=cache->R_ays[j];
  Matrix  R_az=cache->R_azs[j];
  const Matrix & K=cache->Ks[j];
  const double & fx=K.get(0,0);
  const double & fy=K.get(1,1);
  const double * Rd=R.data;
  double scf=1.0;
  if(cache->dofocal)scf=aj[6];
  

  double X=bi[0];
  double Y=bi[1];
  double Z=bi[2];
  double imx= Rd[0]*bi[0]+Rd[1]*bi[1]+Rd[2]*bi[2]+aj[3];
  double imy= Rd[3]*bi[0]+Rd[4]*bi[1]+Rd[5]*bi[2]+aj[4];
  double imz= Rd[6]*bi[0]+Rd[7]*bi[1]+Rd[8]*bi[2]+aj[5];
  double imzz=imz*imz;

  //Finally, the actual derivatives
  Aij[0]=scf*fx*((R_ax[0]*X+R_ax[1]*Y+R_ax[2]*Z)/imz
		 -(R_ax[6]*X+R_ax[7]*Y+R_ax[8]*Z)*(imx/imzz));
  Aij[1]=scf*fx*((R_ay[0]*X+R_ay[1]*Y+R_ay[2]*Z)/imz
		 -(R_ay[6]*X+R_ay[7]*Y+R_ay[8]*Z)*(imx/imzz));
  Aij[2]=scf*fx*((R_az[0]*X+R_az[1]*Y+R_az[2]*Z)/imz
		 -(R_az[6]*X+R_az[7]*Y+R_az[8]*Z)*(imx/imzz));
  Aij[3]=scf*fx*1.0/imz;
  Aij[4]=0;
  Aij[5]=-scf*fx*imx/imzz;
  if(cache->dofocal)Aij[6]=fx*imx/imz;

  int yind=6+cache->dofocal;
  Aij[yind]=scf*fy*((R_ax[3]*X+R_ax[4]*Y+R_ax[5]*Z)/imz
		 -(R_ax[6]*X+R_ax[7]*Y+R_ax[8]*Z)*(imy/imzz));
  Aij[yind+1]=scf*fy*((R_ay[3]*X+R_ay[4]*Y+R_ay[5]*Z)/imz
		 -(R_ay[6]*X+R_ay[7]*Y+R_ay[8]*Z)*(imy/imzz));
  Aij[yind+2]=scf*fy*((R_az[3]*X+R_az[4]*Y+R_az[5]*Z)/imz
		 -(R_az[6]*X+R_az[7]*Y+R_az[8]*Z)*(imy/imzz));
  Aij[yind+3]=0;
  Aij[yind+4]=scf*fy*1.0/imz;
  Aij[yind+5]=-scf*fy*imy/imzz;
  if(cache->dofocal)Aij[yind+6]=fy*imy/imz;
  
  //same as projective case
  Bij[0]=scf*fx*(Rd[0]/imz-Rd[6]*imx/imzz);
  Bij[1]=scf*fx*(Rd[1]/imz-Rd[7]*imx/imzz);
  Bij[2]=scf*fx*(Rd[2]/imz-Rd[8]*imx/imzz);

  Bij[3]=scf*fy*(Rd[3]/imz-Rd[6]*imy/imzz);
  Bij[4]=scf*fy*(Rd[4]/imz-Rd[7]*imy/imzz);
  Bij[5]=scf*fy*(Rd[5]/imz-Rd[8]*imy/imzz);
}

void proj(int j, int i,
	  double *aj, double *bi,
	  double *xij,void *adata){
  double x= aj[0]*bi[0]+aj[1]*bi[1]+aj[2]*bi[2] +aj[3];
  double y= aj[4]*bi[0]+aj[5]*bi[1]+aj[6]*bi[2] +aj[7];
  double z= aj[8]*bi[0]+aj[9]*bi[1]+aj[10]*bi[2]+aj[11];
  xij[0]=x/z;
  xij[1]=y/z;
}
void projac(int j, int i,
	    double *aj, double *bi,
	    double *Aij, double *Bij, void *adata){
  double X=bi[0];
  double Y=bi[1];
  double Z=bi[2];
  double imx= aj[0]*bi[0]+aj[1]*bi[1]+aj[2]*bi[2] +aj[3];
  double imy= aj[4]*bi[0]+aj[5]*bi[1]+aj[6]*bi[2] +aj[7];
  double imz= aj[8]*bi[0]+aj[9]*bi[1]+aj[10]*bi[2]+aj[11];
  double imzz=imz*imz;


  Aij[  0]=X/imz;
  Aij[  1]=Y/imz;
  Aij[  2]=Z/imz;
  Aij[  3]=1.0/imz;
  Aij[  4]=0;
  Aij[  5]=0;
  Aij[  6]=0;
  Aij[  7]=0;
  Aij[  8]=  -X*imx/imzz;
  Aij[  9]=  -Y*imx/imzz;
  Aij[ 10]=  -Z*imx/imzz;
  Aij[ 11]=-1.0*imx/imzz;
  
  Aij[12]=0;
  Aij[13]=0;
  Aij[14]=0;
  Aij[15]=0;
  Aij[16]=X/imz;
  Aij[17]=Y/imz;
  Aij[18]=Z/imz;
  Aij[19]=1.0/imz;
  Aij[20]=  -X*imy/imzz;
  Aij[21]=  -Y*imy/imzz;
  Aij[22]=  -Z*imy/imzz;
  Aij[23]=-1.0*imy/imzz;
  
  Bij[0]=aj[0]/imz-aj[8]*imx/imzz;
  Bij[1]=aj[1]/imz-aj[9]*imx/imzz;
  Bij[2]=aj[2]/imz-aj[10]*imx/imzz;

  Bij[3]=aj[4]/imz-aj[8]*imy/imzz;
  Bij[4]=aj[5]/imz-aj[9]*imy/imzz;
  Bij[5]=aj[6]/imz-aj[10]*imy/imzz;

  //for(int i=0;i<24;i++)Aij[i]*=sign*0;
  //for(int i=0;i<6;i++)Bij[i]*=0;
}

double bundleAdjustSparse(vector<Matrix> & Ps,
			  vector<Matrix> & xs,
			  Matrix & X){
  int n=X.n; //number of points
  int m=Ps.size();//number of images

  char * vis=new char[m*n];//
  double * x=new double[m*n*2];
  
  memset(vis,0,sizeof(char)*m*n);
  memset(x,0,sizeof(double)*m*n*2);

  for(int i=0;i<m;i++){
    for(int j=0;j<X.n;j++){
      vis[j*m+i]=1;

      x[j*m*2+i*2  ]=xs[i](0,j);
      x[j*m*2+i*2+1]=xs[i](1,j);
    }
  }
  int mcon=1;//don't optimize first camera
  double opts[SBA_OPTSSZ]={SBA_INIT_MU,1e-10,1e-10,1e-10};
  double info[SBA_INFOSZ];

  double * adata=0;
  int maxits=100;
  int verbose=0;

  double * p=new double[12*m+3*n];

  //Pack in the data
  for(int i=0;i<m;i++)memcpy(p+12*i,Ps[i].data,sizeof(double)*12);
  
  for(int i=0;i<n;i++){
    p[12*m+3*i  ]=X(0,i)/X(3,i);
    p[12*m+3*i+1]=X(1,i)/X(3,i);
    p[12*m+3*i+2]=X(2,i)/X(3,i);
  }

  int ret=sba_motstr_levmar(n,m,mcon,vis,p,12,3,x,2,
			    proj,projac,
			    adata,maxits,verbose,opts,info);
  if(ret<0)printf(" sba returned %d (error: %f reduced to %f)\n",ret,info[0]/(2.0*m*n),info[1]/(2.0*m*n));
  //Unpack the data
  for(int i=0;i<m;i++){
    memcpy(Ps[i].data,p+12*i,sizeof(double)*12);
    Ps[i]*=(1.0/sqrt(Ps[i].dot(Ps[i])));
  } 
  for(int i=0;i<n;i++){
    X(0,i)=p[12*m+3*i  ];
    X(1,i)=p[12*m+3*i+1];
    X(2,i)=p[12*m+3*i+2];
    X(3,i)=1;
  }
  

  delete [] x;
  delete [] vis;
  delete [] p;

  return info[1]/(2*m*n);//return average reprojection the error
}

void bundleAdjustSparse(vector<Matrix> & Ps,
			vector<vector<Feature *> > & features,
			vector<vector<int> > & x2X,
			vector<Matrix> & X)
{
  const int n=X.size(); //number of points
  const int m=Ps.size();//number of images
  int mcon=0;//think this is right
  //visibility, assuming row major vis[i][j] means point i in image j
  char * vis=new char[m*n];//
  double * xs=new double[m*n*2];

  memset(vis,0,sizeof(char)*m*n);
  memset(xs,0,sizeof(double)*m*n*2);
  
  int xsind=0;
  int nvis=0;
  for(int j=0;j<(int)X.size();j++){
    //for each camera
    for(int i=0;i<m;i++){
      int nvisthiscam=0;
      //for each observation, check to find one that observes this point
      for(int k=0; k<(int)x2X[i].size(); k++){
	int xind=x2X[i][k];
	if(xind==j){
	  vis[xind*m+i]=1;
	  //only add it once
	  if(nvisthiscam==0){
	    xs[xsind++]=features[i][k]->point.x;
	    xs[xsind++]=features[i][k]->point.y;
	    nvis++;
	  }
	  nvisthiscam++;
	}
      }
      //as the point is only added once, we can report problems without
      //the minimization going horribly wrong
      if(nvisthiscam>1){
	fprintf(stderr,"there is a serious problem, a 3D point is visible more than once in the same image (%d times)\n",nvisthiscam);
      }
    }
  }
  if(nvis>m*n){
    printf("this is very very bad!!!!\n");
  }

  double opts[SBA_OPTSSZ]={SBA_INIT_MU,1e-10,1e-10,1e-10};
  double info[SBA_INFOSZ];
  double * adata=0;
  int maxits=100;
  int verbose=0;

  double * p=new double[12*m+3*n];

  //Pack in the data
  for(int i=0;i<m;i++){
    memcpy(p+12*i,Ps[i].data,sizeof(double)*12);
  }
  for(int i=0;i<n;i++){
    p[12*m+3*i  ]=X[i][0]/X[i][3];
    p[12*m+3*i+1]=X[i][1]/X[i][3];
    p[12*m+3*i+2]=X[i][2]/X[i][3];
  }
  // int ret = 
  sba_motstr_levmar(n,m,mcon,vis,p,12,3,xs,2,
                    proj,projac,
                    adata,maxits,verbose,opts,info);
  //printf(" sba returned %d (error: %f reduced to %f)\n",ret,info[0]/(nvis),info[1]/(nvis));
  //Unpack the data
  for(int i=0;i<m;i++){
    memcpy(Ps[i].data,p+12*i,sizeof(double)*12);
    Ps[i]*=(1.0/sqrt(Ps[i].dot(Ps[i])));
  }
  for(int i=0;i<n;i++){
    X[i]=Matrix(4,1);
    X[i][0]=p[12*m+3*i  ];
    X[i][1]=p[12*m+3*i+1];
    X[i][2]=p[12*m+3*i+2];
    X[i][3]=1;
  }
  delete [] xs;
  delete [] vis;
  delete [] p;
}

//euclidean BA
void bundleAdjustSparse(vector<Matrix> & Ks,
			vector<Matrix> & Rts,
			vector<vector<Feature *> > & features,
			vector<vector<int> > & x2X,
			Matrix & X,
			int dofocal)
{
  EuclideanBACache cache(Rts.size(),Ks,dofocal);
  const int n=X.n; //number of points
  const int m=Rts.size();//number of images
  int mcon=0;//think this is right
  //visibility, assuming row major vis[i][j] means point i in image j
  char * vis=new char[m*n];//
  double * xs=new double[m*n*2];

  memset(vis,0,sizeof(char)*m*n);
  memset(xs,0,sizeof(double)*m*n*2);
  
  int xsind=0;
  int nvis=0;
  for(int j=0;j<X.n;j++){
    //for each camera
    for(int i=0;i<m;i++){
      int nvisthiscam=0;
      //for each observation, check to find one that observes this point
      for(int k=0;k<(int)x2X[i].size();k++){
	int xind=x2X[i][k];
	if(xind==j){
	  vis[xind*m+i]=1;
	  //only add it once
	  if(nvisthiscam==0){
	    Matrix obs(3,1);
	    obs[0]=features[i][k]->point.x;
	    obs[1]=features[i][k]->point.y;
	    obs[2]=1;
	    
	    xs[xsind++]=obs[0];
	    xs[xsind++]=obs[1];
	    nvis++;
	  }
	  nvisthiscam++;
	}
      }
      //as the point is only added once, we can report problems without
      //the minimization going horribly wrong
      if(nvisthiscam>1){
	fprintf(stderr,"there is a serious problem in euclidean BA, a 3D point is visible more than once in the same image (%d times)\n",nvisthiscam);
      }
    }
  }
  if(nvis>m*n){
    printf("this is very very bad!!!!\n");
  }

  double opts[SBA_OPTSSZ]={SBA_INIT_MU,1e-10,1e-10,1e-10};
  double info[SBA_INFOSZ];
  void * adata=&cache;
  int maxits=100;
  int verbose=0;

  int cpi=6+dofocal;
  double * p=new double[cpi*m+3*n];
  

  if(2*nvis < cpi*m+3*n){
    fprintf(stderr,"Error: not enough measurements (%d) for problem (%d vars) in %s:%d!",2*nvis,cpi*m+3*n,__FILE__,__LINE__);
    fprintf(stderr,"Error: sba is going to crash !!!!\n");
    
    delete [] p;
    delete [] xs;
    delete [] vis;
    return;
  }


  //Pack in the data
  for(int i=0;i<m;i++){
    Matrix R=Rts[i].submatrix(0,0,3,3);
    Matrix t=Rts[i].submatrix(0,3,3,1);
    Matrix rod=R.matrix2rod();
    memcpy(p+cpi*i,rod.data,sizeof(double)*3);
    memcpy(p+cpi*i+3,t.data,sizeof(double)*3);
    if(dofocal)p[cpi*i+6]=1.0;
  }
  for(int i=0;i<n;i++){
    p[cpi*m+3*i  ]=X(0,i)/X(3,i);
    p[cpi*m+3*i+1]=X(1,i)/X(3,i);
    p[cpi*m+3*i+2]=X(2,i)/X(3,i);
  }
  int ret=sba_motstr_levmar(n,m,mcon,vis,p,cpi,3,xs,2,
			    eucproj,eucprojac,
			    adata,maxits,verbose,opts,info);
  printf(" sba euc returned %d (error: %f reduced to %f)\n",
	 ret,info[0]/(nvis),info[1]/(nvis));
  /*
  // Test with numerical derivatives 
  sba_motstr_levmar(n,m,mcon,vis,p,cpi,3,xs,2,
		    eucproj,0,
		    adata,maxits,verbose,opts,info);
  printf(" sba euc returned %d (error: %f reduced to %f)\n",
  ret,info[0]/(nvis),info[1]/(nvis));*/
  //Unpack the data
  double fx=Ks[0](0,0);
  double fy=Ks[0](1,1);
  double scale_avg=0;
  
  for(int i=0;i<m;i++){
    Matrix rod(3,1);
    Matrix t(3,1);
    memcpy(rod.data,p+cpi*i,sizeof(double)*3);
    memcpy(t.data,p+cpi*i+3,sizeof(double)*3);
    Rts[i].setSubmatrix(0,0,rod.rod2matrix());
    Rts[i].setSubmatrix(0,3,t);
    if(dofocal){
      if(Rts.size()==Ks.size()){
	printf("setting the focal length to %f,%f\n",p[cpi*i+6]*fx,p[cpi*i+6]*fy);
	Ks[i](0,0)=p[cpi*i+6]*fx;
	Ks[i](1,1)=p[cpi*i+6]*fy;	
      }
      else{
	printf("FIXME: focal %d is %f of (%f,%f)=>(%f,%f)\n",
	       i,p[cpi*i+6],fx,fy,p[cpi*i+6]*fx,p[cpi*i+6]*fy);
	scale_avg+=p[cpi*i+6];
      }
    }
  }
  if(dofocal && Rts.size()!=Ks.size()){
    scale_avg/=m;
    printf("There are not enough intrinsics, updating the %ld existing to be the same (scale_avg=%f)\n",
           Ks.size(), scale_avg);
    for(int i=0; i<(int)Ks.size();i++){
      Ks[i](0,0)=fx*scale_avg;
      Ks[i](1,1)=fy*scale_avg;
    }
  }
  
  for(int i=0;i<n;i++){
    X(0,i)=p[cpi*m+3*i  ];
    X(1,i)=p[cpi*m+3*i+1];
    X(2,i)=p[cpi*m+3*i+2];
    X(3,i)=1;
  }
  delete [] xs;
  delete [] vis;
  delete [] p;
}

namespace BAEuclideanCamera{
  Matrix xs;
  Matrix Xg;
  Matrix K;
  double g_err;
  void proj(const int & m,const int & n,
	    double * aj,double * fvec,int & info){
    Matrix rod(3,1);
    rod[0]=aj[0];
    rod[1]=aj[1];
    rod[2]=aj[2];
    Matrix R=rod.rod2matrix();
    const double * Rd=R.data;
    double fx=K.get(0,0);
    double fy=K.get(1,1);
    double px=K.get(0,2);
    double py=K.get(1,2);
    double scf=1.0;
    //int dofocal=(n==7);
    if(n==7)scf=aj[6];
    
    g_err=0;
    for(int i=0;i<Xg.n;i++){
      double X=Xg(0,i);
      double Y=Xg(1,i);
      double Z=Xg(2,i);
      double imx= (Rd[0]*X+Rd[1]*Y+Rd[2]*Z+aj[3]);
      double imy= (Rd[3]*X+Rd[4]*Y+Rd[5]*Z+aj[4]);
      double imz= (Rd[6]*X+Rd[7]*Y+Rd[8]*Z+aj[5]);
      fvec[2*i  ]=(scf*fx*imx/imz+px)-xs(0,i);
      fvec[2*i+1]=(scf*fy*imy/imz+py)-xs(1,i);
      g_err+=fvec[2*i]*fvec[2*i];
      g_err+=fvec[2*i+1]*fvec[2*i+1];;
    }
    g_err/=(m);
  }
  void projac(const int & m,const int & n,
	      double * aj,double * fvec,double * fjac,
	      const int & ldfjac,int & info){
     Matrix rod(3,1);
     rod[0]=aj[0];
     rod[1]=aj[1];
     rod[2]=aj[2];
     Matrix R=rod.rod2matrix();
     const double * Rd=R.data;
     double fx=K.get(0,0);
     double fy=K.get(1,1);
     double px=K.get(0,2);
     double py=K.get(1,2);
     double scf=1.0;
     int dofocal=(n==7);
     if(n==7)scf=aj[6];
     if(info==1){
       g_err=0;
       for(int i=0;i<Xg.n;i++){
	 double X=Xg(0,i);
	 double Y=Xg(1,i);
	 double Z=Xg(2,i);
	 double imx= (Rd[0]*X+Rd[1]*Y+Rd[2]*Z+aj[3]);
	 double imy= (Rd[3]*X+Rd[4]*Y+Rd[5]*Z+aj[4]);
	 double imz= (Rd[6]*X+Rd[7]*Y+Rd[8]*Z+aj[5]);
	 fvec[2*i  ]=(scf*fx*imx/imz+px)-xs(0,i);
	 fvec[2*i+1]=(scf*fy*imy/imz+py)-xs(1,i);
	 g_err+=fvec[2*i]*fvec[2*i];
	 g_err+=fvec[2*i+1]*fvec[2*i+1];;
       }
       g_err/=(m);
     }
     if(info==2){
       double angle=sqrt(aj[0]*aj[0]+aj[1]*aj[1]+aj[2]*aj[2]);
       Matrix S(3,3);
       S(0,0)=0;
       S(0,1)=-aj[2]/angle;
       S(0,2)= aj[1]/angle;
       S(1,0)= aj[2]/angle;
       S(1,1)= 0;
       S(1,2)=-aj[0]/angle;
       S(2,0)=-aj[1]/angle;
       S(2,1)= aj[0]/angle;
       S(2,2)=0;
       Matrix SS=S*S;
       
       double sin_angle=sin(angle);
       double cos_angle=cos(angle);
       double angle_ax=aj[0]/angle;
       double angle_ay=aj[1]/angle;
       double angle_az=aj[2]/angle;
       Matrix S_ax(3,3),S_ay(3,3),S_az(3,3);
       Matrix R_ax,R_ay,R_az;
      
       S_ax.setAll(0);
       S_ay.setAll(0);
       S_az.setAll(0);
       S_ax(1,2)=-1;
       S_ax(2,1)= 1;
       S_ay(0,2)=1;
       S_ay(2,0)=-1;
       S_az(0,1)=-1;
       S_az(1,0)=1;
       S_ax=S_ax*(1.0/angle)-S*(angle_ax/angle);
       S_ay=S_ay*(1.0/angle)-S*(angle_ay/angle);
       S_az=S_az*(1.0/angle)-S*(angle_az/angle);
       
       R_ax=S_ax*sin_angle+S*(cos_angle*angle_ax)+
	 (S_ax*S+S*S_ax)*(1.0-cos_angle)+SS*sin_angle*angle_ax;
       R_ay=S_ay*sin_angle+S*(cos_angle*angle_ay)+
	 (S_ay*S+S*S_ay)*(1.0-cos_angle)+SS*sin_angle*angle_ay;
       R_az=S_az*sin_angle+S*(cos_angle*angle_az)+
	 (S_az*S+S*S_az)*(1.0-cos_angle)+SS*sin_angle*angle_az;
       
       for(int i=0;i<Xg.n;i++){
	 double X=Xg(0,i);
	 double Y=Xg(1,i);
	 double Z=Xg(2,i);
	 double imx= Rd[0]*X+Rd[1]*Y+Rd[2]*Z+aj[3];
	 double imy= Rd[3]*X+Rd[4]*Y+Rd[5]*Z+aj[4];
	 double imz= Rd[6]*X+Rd[7]*Y+Rd[8]*Z+aj[5];
	 double imzz=imz*imz;
	 

	 //Finally, the actual derivatives
	 fjac[2*i         ]=scf*fx*((R_ax[0]*X+R_ax[1]*Y+R_ax[2]*Z)/imz
				   -(R_ax[6]*X+R_ax[7]*Y+R_ax[8]*Z)*(imx/imzz));
	 fjac[2*i+ldfjac  ]=scf*fx*((R_ay[0]*X+R_ay[1]*Y+R_ay[2]*Z)/imz
				   -(R_ay[6]*X+R_ay[7]*Y+R_ay[8]*Z)*(imx/imzz));
	 fjac[2*i+2*ldfjac]=scf*fx*((R_az[0]*X+R_az[1]*Y+R_az[2]*Z)/imz
				   -(R_az[6]*X+R_az[7]*Y+R_az[8]*Z)*(imx/imzz));
	 fjac[2*i+3*ldfjac]=scf*fx*(1.0/imz);
	 fjac[2*i+4*ldfjac]=0;
	 fjac[2*i+5*ldfjac]=-scf*fx*(imx/imzz);
	 if(dofocal)fjac[2*i+6*ldfjac]=fx*(imx/imz);

	 fjac[2*i+1         ]=scf*fy*((R_ax[3]*X+R_ax[4]*Y+R_ax[5]*Z)/imz
				      -(R_ax[6]*X+R_ax[7]*Y+R_ax[8]*Z)*(imy/imzz));
	 fjac[2*i+1+ldfjac  ]=scf*fy*((R_ay[3]*X+R_ay[4]*Y+R_ay[5]*Z)/imz
				      -(R_ay[6]*X+R_ay[7]*Y+R_ay[8]*Z)*(imy/imzz));
	 fjac[2*i+1+2*ldfjac]=scf*fy*((R_az[3]*X+R_az[4]*Y+R_az[5]*Z)/imz
				      -(R_az[6]*X+R_az[7]*Y+R_az[8]*Z)*(imy/imzz));
	 fjac[2*i+1+3*ldfjac]=0;
	 fjac[2*i+1+4*ldfjac]=scf*fy*(1.0/imz);
	 fjac[2*i+1+5*ldfjac]=-scf*fy*(imy/imzz);
	 if(dofocal)fjac[2*i+1+6*ldfjac]=fy*(imy/imz);
       }
     }
  }
}

//no need to be sparse for one camera (NOTE: hasn't been tested yet
double bundleAdjustCamera(const Matrix & K,
			  Matrix & Rt,
			  const Matrix & xs,
			  const Matrix & X,int dofocal){
  BAEuclideanCamera::xs=xs.copy();
  BAEuclideanCamera::Xg=Matrix(3,X.n);
  BAEuclideanCamera::K=K;  
  for(int i=0;i<X.n;i++){
    double denom=1;
    if(X.m>3)denom=X.get(3,i);
    BAEuclideanCamera::Xg(0,i)=X.get(0,i)/denom;
    BAEuclideanCamera::Xg(1,i)=X.get(1,i)/denom;
    BAEuclideanCamera::Xg(2,i)=X.get(2,i)/denom;
  }
  
  Matrix x(6+dofocal,1);
  Matrix fvec(2*X.n,1);

  //pack initial solution
  Matrix rod=Rt.submatrix(0,0,3,3).matrix2rod();
  rod.printMatlab("rod");
  x[0]=rod[0];
  x[1]=rod[1];
  x[2]=rod[2];
  x[3]=Rt(0,3);
  x[4]=Rt(1,3);
  x[5]=Rt(2,3);
  if(dofocal)x[6]=1;//scale of existing focal length

  //April 19: residuals were 1e-10 before, way too big
  Matrix::lsqnonlin(BAEuclideanCamera::projac,x,fvec,1e-8,1e-8,1e-8,20000);

  //unpack results;
  rod[0]=x[0];
  rod[1]=x[1];
  rod[2]=x[2];
  Rt.setSubmatrix(0,0,rod.rod2matrix());
  Rt(0,3)=x[3];
  Rt(1,3)=x[4];
  Rt(2,3)=x[5];
  if(dofocal)return x[6];
  return 1.0;
}

//no need to be sparse for one camera
double bundleAdjustCamera(const Matrix & K,
			Matrix & Rt,
			vector<Feature *> & features,
			vector<int> & x2X,
			const Matrix & X,int dofocal){
  vector<int> obs_ind;
  char * xcnt=new char[X.n];
  memset(xcnt,0,sizeof(char)*X.n);
  for(int i=0;i<(int)features.size();i++){
    int xind=x2X[i];
    //only add the first one
    if(xind>=0 && !xcnt[xind]){
      xcnt[xind]=1;
      obs_ind.push_back(i);
    }
  }
  delete [] xcnt;

  BAEuclideanCamera::xs=Matrix(2,obs_ind.size());
  BAEuclideanCamera::Xg=Matrix(3,obs_ind.size());
  BAEuclideanCamera::K=K;  
  for(int i=0;i<(int)obs_ind.size();i++){
    BAEuclideanCamera::xs(0,i)=features[obs_ind[i]]->point.x;
    BAEuclideanCamera::xs(1,i)=features[obs_ind[i]]->point.y;

    int xind=x2X[obs_ind[i]];
    BAEuclideanCamera::Xg(0,i)=X.get(0,xind)/X.get(3,xind);
    BAEuclideanCamera::Xg(1,i)=X.get(1,xind)/X.get(3,xind);
    BAEuclideanCamera::Xg(2,i)=X.get(2,xind)/X.get(3,xind);
  }
    
  Matrix x(6+dofocal,1);
  Matrix fvec(2*obs_ind.size(),1);

  //pack initial solution
  Matrix rod=Rt.submatrix(0,0,3,3).matrix2rod();
  x[0]=rod[0];
  x[1]=rod[1];
  x[2]=rod[2];
  x[3]=Rt(0,3);
  x[4]=Rt(1,3);
  x[5]=Rt(2,3);
  if(dofocal)x[6]=1;//scale of existing focal length
  
  //BAEuclideanCamera::projac(fvec.m,x.m,x.data,fvec.data,0,0,info);
  //printf("BACam before :%f\n",BAEuclideanCamera::g_err);
  Matrix::lsqnonlin(BAEuclideanCamera::projac,x,fvec,1e-8,1e-8,1e-8);

  //BAEuclideanCamera::projac(fvec.m,x.m,x.data,fvec.data,0,0,info);
  //printf("BACam after  :%f\n",BAEuclideanCamera::g_err);

  //unpack results;
  rod[0]=x[0];
  rod[1]=x[1];
  rod[2]=x[2];
  Rt.setSubmatrix(0,0,rod.rod2matrix());
  Rt(0,3)=x[3];
  Rt(1,3)=x[4];
  Rt(2,3)=x[5];
  if(dofocal)return x[6];
  return 1.0;
}
void bundleAdjust(vector<Matrix> & Ps,
		  vector<vector<Feature *> > & features,
		  vector<vector<int> > & x2X,
		  vector<Matrix> & X)
{
  int nobs=0;
  
  BA::npts=X.size();
  BA::ncams=Ps.size();
  BA::xs=vector<Matrix>();
  BA::obs=vector<vector<int> >();

  for(int i=0;i<(int)Ps.size();i++){
    vector<int> find;
    vector<int> ob;
    for(int j=0;j<(int)x2X[i].size();j++){
      if(x2X[i][j]>=0){
	ob.push_back(x2X[i][j]);
	find.push_back(j);
      }
    }
    Matrix x(2,find.size());
    for(int j=0;j<(int)find.size();j++){
      x(0,j)=features[i][find[j]]->point.x;
      x(1,j)=features[i][find[j]]->point.y;
    }
    BA::xs.push_back(x);
    BA::obs.push_back(ob);
    nobs+=ob.size();
  }

  Matrix fvec(2*nobs,1);
  Matrix xvec(12*Ps.size()+3*X.size(),1);
  
  for(int i=0;i<(int)Ps.size();i++)
    memcpy(xvec.data+12*i,Ps[i].data,sizeof(double)*12);
  for(int i=0;i<(int)X.size();i++){
    xvec[12*Ps.size()+3*i  ]=X[i][0]/X[i][3];
    xvec[12*Ps.size()+3*i+1]=X[i][1]/X[i][3];
    xvec[12*Ps.size()+3*i+2]=X[i][2]/X[i][3];    
  }
  int info=1;
  BA::bundleAdjustVis(fvec.m,xvec.m,xvec.data,fvec.data,0,0,info);
  printf("error before %f (ncams:%d,npts:%d,nobs:%d,m:%d,n:%d)\n",BA::g_err,
	 BA::ncams,BA::npts,nobs,fvec.m,xvec.m);

  Matrix::lsqnonlin(BA::bundleAdjustVis,xvec,fvec,
                    1e-6,1e-6,1e-6,20);
  Matrix Xmat;
  BA::unpack(xvec.data,Ps,Xmat);
  for(int j=0;j<(int)X.size();j++)
    X[j]=Xmat.getColumn(j);

  BA::bundleAdjustVis(fvec.m,xvec.m,xvec.data,fvec.data,0,0,info);
  printf("      after %f\n",BA::g_err);
}
		  
double bundleAdjust(vector<Matrix> & Ps,
		    vector<Matrix> & xs,
		    Matrix & X)
{
  if(0)
  {
    printf("Testing the vis BA for quality and speed\n");
    
    vector<vector<Feature *> > features;
    vector<vector<int> > x2X;
    vector<Matrix> Xvec;
    vector<Matrix> Ps2;
    for(int i=0;i<(int)Ps.size();i++){
      Ps2.push_back(Ps[i].copy());
      vector<int> obs;
      vector<Feature * > feats;
      for(int j=0;j<X.n;j++){
	obs.push_back(j);
	feats.push_back(new Feature(0,Vec2f(xs[i](0,j),xs[i](1,j))));
      }
      features.push_back(feats);
      x2X.push_back(obs);
    }
    for(int i=0;i<X.n;i++)
      Xvec.push_back(X.getColumn(i));

    bundleAdjust(Ps2,features,x2X,Xvec);

    printf("done testing the vis BA for quality and speed\n");
    for(int i=0;i<(int)features.size();i++){
      for(int j=0;j<(int)features[i].size();j++)
	delete features[i][j];
    }
  }
  BA::xs=xs;
  BA::ncams=Ps.size();
  BA::npts=X.n;
  printf("doing bundle adjustment on %d cameras with %d npts\n",BA::ncams,BA::npts);

  Matrix x(12*BA::ncams+3*BA::npts,1);
  Matrix fvec(BA::ncams*BA::npts*2,1);

  BA::pack(x.data,Ps,X);

  int info=1;
  BA::bundleAdjustErr(fvec.m,x.m,x.data,fvec.data,info);
  printf("error before %f (m:%d,n:%d)\n",BA::g_err,fvec.m,x.m);

  for(int k=0;k<1;k++){
    Matrix::lsqnonlin(BA::bundleAdjustErrAn,x,fvec,
                      1e-6,1e-6,1e-6,20);
    BA::unpack(x.data,Ps,X);
    //printf("return %d\n",ret);    
  }
  BA::bundleAdjustErr(fvec.m,x.m,x.data,fvec.data,info);
  printf("error after analytical %f\n",BA::g_err);

  return BA::g_err;//return error
  /*
  Matrix::lsqnonlin(BA::bundleAdjustErr,x,fvec);
  BA::unpack(x.data,Ps,X);

  BA::bundleAdjustErr(fvec.m,x.m,x.data,fvec.data,info);
  printf("error after numerical %f\n",BA::g_err);*/
}
