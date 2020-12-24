#include "recon_globals.h"
#include "recon_geometry.h"
#include "bundle_adjust.h"
#include  <sys/time.h>
#include <algorithm>
#include <assert.h>
#include <nmisc/timer.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>

using namespace std;
using nacb::Vec3;
using nacb::Vec2;
using nacb::Vec3d;

Matrix rectifyPoints(Matrix & p)
{
  Matrix B(3,3);
  Matrix A(3,3);
  for(int i=0;i<3;i++){
    A(0,i)=p(0,i);
    A(1,i)=p(1,i);
    A(2,i)=1;
  }
  Matrix x3(3,1);
  x3[0]=p(0,3);
  x3[1]=p(1,3);
  x3[2]=1;
  
  Matrix lam=A.inverse()*x3;

  for(int i=0;i<3;i++){
    //double len=sqrt(p(0,i)*p(0,i)+p(1,i)*p(1,i));
    //lam[i]=1;
    B(0,i)=lam[i]*p(0,i);
    B(1,i)=lam[i]*p(1,i);
    B(2,i)=lam[i];
  }
  B=B.inverse();
  return B;
}


template <class T>
Vec3<T> triangulate(const Matrix & P1, const Matrix & P2,
		    const Vec2<T> & p1, const Vec2<T> & p2){
  Matrix A(4,3);

  A(0,0)=P1(0,0)-P1(2,0)*p1.x;
  A(0,1)=P1(0,1)-P1(2,1)*p1.x;
  A(0,2)=P1(0,2)-P1(2,2)*p1.x;
  A(1,0)=P1(1,0)-P1(2,0)*p1.y;
  A(1,1)=P1(1,1)-P1(2,1)*p1.y;
  A(1,2)=P1(1,2)-P1(2,2)*p1.y;
  A(2,0)=P2(0,0)-P2(2,0)*p2.x;
  A(2,1)=P2(0,1)-P2(2,1)*p2.x;
  A(2,2)=P2(0,2)-P2(2,2)*p2.x;
  A(3,0)=P2(1,0)-P2(2,0)*p2.y;
  A(3,1)=P2(1,1)-P2(2,1)*p2.y;
  A(3,2)=P2(1,2)-P2(2,2)*p2.y;
  
  Matrix b(4,1);
  b[0]=P1(2,3)*p1.x-P1(0,3);
  b[1]=P1(2,3)*p1.y-P1(1,3);
  b[2]=P2(2,3)*p2.x-P2(0,3);
  b[3]=P2(2,3)*p2.y-P2(1,3);
  
  Matrix X=Matrix::LlinLeastSq(A,b);
  return Vec3<T>(X[0],X[1],X[2]);
}


Matrix triangulate(vector<Matrix> & Ps,
		   vector<Matrix> & xs)
{
  assert(Ps.size()>=2);
  Matrix Xs(4,xs[0].n);
  
  for(int i=0;i<(int)Ps.size();i++){
    for(int j=0;j<12;j++){
      if(isinf(Ps[i][j])||
	 isnan(Ps[i][j])){
	printf("it is nan, in triangulate\n");
	return Xs;
      }
    }
  }
  
  //for each point
  for(int j=0;j<Xs.n;j++){
    //for each camera
    Matrix A(Ps.size()*2,4);
    for(int i=0;i<(int)Ps.size();i++){
      for(int k=0;k<4;k++){
	A(2*i,  k)=xs[i](0,j)*Ps[i](2,k)-Ps[i](0,k);
	A(2*i+1,k)=xs[i](1,j)*Ps[i](2,k)-Ps[i](1,k);
      }
    }
    Matrix U,S,V;
    A.Lsvd(U,S,V);

    // S.printMatlab("S");

    Matrix X=V.getColumn(3);
    Xs(0,j)=X[0]/X[3];
    Xs(1,j)=X[1]/X[3];
    Xs(2,j)=X[2]/X[3];
    Xs(3,j)=1;
    //(A*X).printMatlab("A*X");
  }
  return Xs;
}
//hasn't been tested yet
Matrix resect(Matrix & X,Matrix & x){
  Matrix A(X.n*2,12);
  Matrix P=Matrix::eye(3,4);
  for(int its=0;its<10;its++){
    for(int i=0;i<X.n;i++){
      double w=1;
      if(its){
	Matrix pp=P*X.getColumn(i);
	w=1.0/fabs(pp[2]+1e-10);
      }
      A(2*i+0,0)=X(0,i)*w;
      A(2*i+0,1)=X(1,i)*w;
      A(2*i+0,2)=X(2,i)*w;
      A(2*i+0,3)=X(3,i)*w;
      A(2*i+0,4)=0;
      A(2*i+0,5)=0;
      A(2*i+0,6)=0;
      A(2*i+0,7)=0;
      A(2*i+0, 8)=-X(0,i)*x(0,i)*w;
      A(2*i+0, 9)=-X(1,i)*x(0,i)*w;
      A(2*i+0,10)=-X(2,i)*x(0,i)*w;
      A(2*i+0,11)=-X(3,i)*x(0,i)*w;
    
      A(2*i+1,0)=0;
      A(2*i+1,1)=0;
      A(2*i+1,2)=0;
      A(2*i+1,3)=0;
      A(2*i+1,4)=X(0,i)*w;
      A(2*i+1,5)=X(1,i)*w;
      A(2*i+1,6)=X(2,i)*w;
      A(2*i+1,7)=X(3,i)*w;
      A(2*i+1, 8)=-X(0,i)*x(1,i)*w;
      A(2*i+1, 9)=-X(1,i)*x(1,i)*w;
      A(2*i+1,10)=-X(2,i)*x(1,i)*w;
      A(2*i+1,11)=-X(3,i)*x(1,i)*w;
    }
    Matrix U,S,V;
    if(A.m >= 1000)
      (A.transpose()*A).Lsvd(U,S,V);
    else
      A.Lsvd(U,S,V);
    memcpy(P.data,(V.getColumn(V.n-1)).data,sizeof(double)*12);
  }
  return P;
}
//hasn't been tested yet
Matrix resectEuclidean(Matrix & K,Matrix & X,Matrix & x,bool dofocal){
  Matrix P=resect(X,x);
  Matrix Rt=K.inverse()*P;
  
  Matrix U,S,V;
  (Rt.submatrix(0,0,3,3)).Lsvd(U,S,V);
  Rt.setSubmatrix(0,0,U*V.transpose());

  /*double fsc=*/
  bundleAdjustCamera(K,Rt,x,X,dofocal);
  
  return Rt;
}

//hasn't been tested yet
vector<int> resectEuclideanRobust(Matrix & Rt,Matrix & K,Matrix & X,Matrix & x){
  int cset[6];
  for(int its=0;its<100;its++){
    unique_set(cset,X.n,6);
    Matrix xuse(2,6);
    Matrix Xuse(4,6);
    for(int i=0;i<6;i++){
      xuse(0,i)=x(0,cset[i]);
      xuse(1,i)=x(1,cset[i]);

      Xuse(0,i)=X(0,cset[i]);
      Xuse(1,i)=X(1,cset[i]);
      Xuse(2,i)=X(2,cset[i]);
      Xuse(3,i)=X(3,cset[i]);
    }
    
  }
  printf("this isn't done yet!!!!!\n");
  return vector<int>();
}

namespace Align{
  Matrix xs;
  vector<int>   cam;
  vector<Matrix> Ps;
  Matrix X;
  double g_err;
  void alignErr(const int & m,const int & n,
		double * x,double * fvec,
		double * fjac,const int & ldfjac,
		int & info){
    Matrix H(4,4);
    Matrix HX;
    memcpy(H.data,x,sizeof(double)*16);
    HX=H*X;
    
    vector<Matrix> phxs;
    vector<Matrix> PHXs;
    for(int i=0;i<(int)Ps.size();i++){
      PHXs.push_back(Ps[i]*HX);
      phxs.push_back(projectPoints(Ps[i],HX));
    }

    if(info==1){
      double err=0;
      for(int j=0;j<xs.n;j++){
	fvec[2*j]=  phxs[cam[j]](0,j)-xs(0,j);
	fvec[2*j+1]=phxs[cam[j]](1,j)-xs(1,j);
	err+=fvec[2*j]*fvec[2*j];
	err+=fvec[2*j+1]*fvec[2*j+1];
      }
      if(isnan(err)||isinf(err))info=-1;
      //printf("Align: %f\n",err);
    }
    if(info==2){
      for(int j=0;j<xs.n;j++){
	int xind=j;
	
	Matrix & P=Ps[cam[j]];
	Matrix & PHX=PHXs[cam[j]];

	for(int ii=0;ii<4;ii++){
	  for(int jj=0;jj<4;jj++){
	    fjac[(2*j  )+(4*ii+jj)*ldfjac]=P(0,ii)*X(jj,xind)/PHX(2,xind)-
	      (PHX(0,xind)/(PHX(2,xind)*PHX(2,xind)))*P(2,ii)*X(jj,xind);
	    fjac[(2*j+1)+(4*ii+jj)*ldfjac]=P(1,ii)*X(jj,xind)/PHX(2,xind)-
	      (PHX(1,xind)/(PHX(2,xind)*PHX(2,xind)))*P(2,ii)*X(jj,xind);
	  }
	}
      }
    }
  }
  void alignEucErrFD(const int & m,const int & n,
		     double * x,double * fvec,
		     int & info){
    Matrix rod(3,1);
    rod[0]=x[0];
    rod[1]=x[1];
    rod[2]=x[2];
    Matrix R=rod.rod2matrix();
    Matrix H=Matrix::eye(4,4);
    H.setSubmatrix(0,0,R);
    H(0,3)=x[3];
    H(1,3)=x[4];
    H(2,3)=x[5];

    for(int i=0;i<3;i++){
      for(int j=0;j<4;j++)
	H(i,j)*=x[6];
    }
    Matrix HX=H*X;
    vector<Matrix> phxs;
    vector<Matrix> PHXs;
    for(int i=0;i<(int)Ps.size();i++){
      PHXs.push_back(Ps[i]*HX);
      phxs.push_back(projectPoints(Ps[i],HX));
    }

    double err=0;
    for(int j=0;j<xs.n;j++){
      fvec[2*j]=  phxs[cam[j]](0,j)-xs(0,j);
      fvec[2*j+1]=phxs[cam[j]](1,j)-xs(1,j);
      err+=fvec[2*j]*fvec[2*j];
      err+=fvec[2*j+1]*fvec[2*j+1];
    }
    if(isnan(err)||isinf(err))info=-1;
    g_err=err/(2*xs.n);
    //printf("Align: %f\n",err);
  }
  void alignEucErr(const int & m,const int & n,
		   double * x,double * fvec,
		   double * fjac,const int & ldfjac,
		   int & info){
    Matrix rod(3,1);
    rod[0]=x[0];
    rod[1]=x[1];
    rod[2]=x[2];
    Matrix R=rod.rod2matrix();
    Matrix H=Matrix::eye(4,4);

    H.setSubmatrix(0,0,R);
    H(0,3)=x[3];
    H(1,3)=x[4];
    H(2,3)=x[5];
    Matrix dHdsc=H.copy();
    dHdsc(3,3)=0;

    for(int i=0;i<3;i++){
      for(int j=0;j<4;j++)
	H(i,j)*=x[6];
    }
    Matrix HX=H*X;
    vector<Matrix> phxs;
    vector<Matrix> PHXs;
    for(int i=0;i<(int)Ps.size();i++){
      PHXs.push_back(Ps[i]*HX);
      phxs.push_back(projectPoints(Ps[i],HX));
    }
    if(info==1){
      double err=0;
      for(int j=0;j<xs.n;j++){
	fvec[2*j]=  phxs[cam[j]](0,j)-xs(0,j);
	fvec[2*j+1]=phxs[cam[j]](1,j)-xs(1,j);
	err+=fvec[2*j]*fvec[2*j];
	err+=fvec[2*j+1]*fvec[2*j+1];
      }
      if(isnan(err)||isinf(err))info=-1;
      g_err=err/(2*xs.n);
      //printf("Align: %f\n",err);
    }
    if(info==2){
      double angle=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      Matrix S(3,3);
      S(0,0)=0;
      S(0,1)=-x[2]/angle;
      S(0,2)= x[1]/angle;
      S(1,0)= x[2]/angle;
      S(1,1)= 0;
      S(1,2)=-x[0]/angle;
      S(2,0)=-x[1]/angle;
      S(2,1)= x[0]/angle;
      S(2,2)=0;
      Matrix SS=S*S;
      
      double sin_angle=sin(angle),cos_angle=cos(angle);
      double angle_ax=x[0]/angle,angle_ay=x[1]/angle,angle_az=x[2]/angle;
      Matrix S_ax(3,3),S_ay(3,3),S_az(3,3);
      Matrix dHdx[7],R_axyz[3];
      
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
      
      R_axyz[0]=S_ax*sin_angle+S*(cos_angle*angle_ax)+
	(S_ax*S+S*S_ax)*(1.0-cos_angle)+SS*sin_angle*angle_ax;
      R_axyz[1]=S_ay*sin_angle+S*(cos_angle*angle_ay)+
	(S_ay*S+S*S_ay)*(1.0-cos_angle)+SS*sin_angle*angle_ay;
      R_axyz[2]=S_az*sin_angle+S*(cos_angle*angle_az)+
	(S_az*S+S*S_az)*(1.0-cos_angle)+SS*sin_angle*angle_az;

      for(int k=0;k<7;k++){
	dHdx[k]=Matrix(4,4);
	dHdx[k].setAll(0);
      }
      dHdx[0].setSubmatrix(0,0,R_axyz[0]*x[6]);
      dHdx[1].setSubmatrix(0,0,R_axyz[1]*x[6]);
      dHdx[2].setSubmatrix(0,0,R_axyz[2]*x[6]);
      
      dHdx[3](0,3)=x[6];
      dHdx[4](1,3)=x[6];
      dHdx[5](2,3)=x[6];

      dHdx[6]=dHdsc;

      vector<vector<Matrix> > Ps_dHdxs;
      for(int i=0;i<(int)Ps.size();i++){
	vector<Matrix> P_dHdxs;
	for(int k=0;k<7;k++){
	  P_dHdxs.push_back(Ps[i]*dHdx[k]);
	}
	Ps_dHdxs.push_back(P_dHdxs);
      }
      for(int j=0;j<xs.n;j++){
	int xind=j;
	
	Matrix & PHX=PHXs[cam[j]];
	Matrix Xxind=X.getColumn(xind);
	
	//this could be faster, because lots of dHdx is 0,
	//but this is soo much easier to read
	for(int k=0;k<7;k++){
	  Matrix P_dHdx_X=Ps_dHdxs[cam[j]][k]*Xxind;
	  fjac[(2*j  )+ldfjac*k]=P_dHdx_X[0]/PHX(2,xind)
	    -PHX(0,xind)/(PHX(2,xind)*PHX(2,xind))*P_dHdx_X[2];
	  fjac[(2*j+1)+ldfjac*k]=P_dHdx_X[1]/PHX(2,xind)
	    -PHX(1,xind)/(PHX(2,xind)*PHX(2,xind))*P_dHdx_X[2];
	}
      }
      /*
	//checking jacobian
      double * fvec2=new double[m];
      
      for(int k=0;k<7;k++){
	int info2=1;
	Matrix x2(7,1);
	memcpy(x2.data,x,sizeof(double)*7);
	double dx=1e-7;
	x2[k]+=dx;
	alignEucErr(m,n,x2.data,fvec2,0,0,info2);
	for(int h=0;h<m;h++){
	  double shouldbe=(fvec2[h]-fvec[h])/dx;
	  printf("%f %f\n",shouldbe,fjac[h+ldfjac*k]);
	}
      }
      delete [] fvec2;*/
    }
  }
};

void align_refine(Matrix & H,vector<Matrix> & Ps,
		  Matrix & X,Matrix & xs,vector<int> & cams)
{
  Align::xs=xs;
  Align::X=X;
  Align::cam=cams;
  Align::Ps=Ps;

  Matrix fvec(2*xs.n,1);
  Matrix x(16,1);

  memcpy(x.data,H.data,sizeof(double)*16);
  
  Matrix::lsqnonlin(Align::alignErr,x,fvec);
  
  memcpy(H.data,x.data,sizeof(double)*16);
}

void align_refine_euc(Matrix & R,Matrix & t,double & sc,
		      vector<Matrix> & Ps,
		      Matrix & X,Matrix & xs,vector<int> & cams)
{
  Align::xs=xs;
  Align::X=X;
  Align::cam=cams;
  Align::Ps=Ps;

  Matrix fvec(2*xs.n,1);
  Matrix x(7,1);

  Matrix rod=R.matrix2rod();
  x[0]=rod[0];
  x[1]=rod[1];
  x[2]=rod[2];
  x[3]=t[0];
  x[4]=t[1];
  x[5]=t[2];
  x[6]=sc;

  int info=1;
  Align::alignEucErr(fvec.m,x.m,x.data,fvec.data,0,0,info);
  //double before=Align::g_err;

  Matrix::lsqnonlin(Align::alignEucErr,x,fvec,
                    1e-5,1e-6,1e-5,1000);

  info=1;Align::alignEucErr(fvec.m,x.m,x.data,fvec.data,0,0,info);
  //printf(" euc align %f --> %f\n",before,Align::g_err);

  rod[0]=x[0];
  rod[1]=x[1];
  rod[2]=x[2];
  R=rod.rod2matrix();
  t[0]=x[3];
  t[1]=x[4];
  t[2]=x[5];
  sc=x[6];
}

void alignImageBased(Matrix & H,vector<Matrix> & Ps,
		     Matrix & X,Matrix & x,vector<int> & cams,
		     int affine=0){
  int nvars=affine?12:16;
  H=Matrix::eye(4,4);
  Matrix A(X.n*2,nvars);

  for(int its=0;its<20;its++){
    vector<Matrix> PHs;
    
    for(int i=0;i<(int)Ps.size();i++){
      Matrix P=Ps[i];
      Matrix PH=P*H;
      PH*=(1.0/sqrt(PH.dot(PH)));
      PHs.push_back(PH);
    }

    for(int i=0;i<X.n;i++){
      Matrix & P=Ps[cams[i]];
      Matrix & PH=PHs[cams[i]];
      Matrix X2=X.getColumn(i);
	  
      double u=x(0,i);
      double v=x(1,i);
      
      int ind1=2*i;
      int ind2=2*i+1;
	  
      double w=1;
      if(its!=0){
	Matrix pph=PH*X2;
	w=1.0/(fabs(pph[2])+1e-10);
      }
      //else{
      //Matrix pph=P*X[x2X[i1][j]];
      //w=1.0/pph[2];
      //}
	  
      A(ind1,0)=(P(0,0)*X2[0]-u*P(2,0)*X2[0])*w;
      A(ind1,1)=(P(0,0)*X2[1]-u*P(2,0)*X2[1])*w;
      A(ind1,2)=(P(0,0)*X2[2]-u*P(2,0)*X2[2])*w;
      A(ind1,3)=(P(0,0)*X2[3]-u*P(2,0)*X2[3])*w;
      A(ind1,4)=(P(0,1)*X2[0]-u*P(2,1)*X2[0])*w;
      A(ind1,5)=(P(0,1)*X2[1]-u*P(2,1)*X2[1])*w;
      A(ind1,6)=(P(0,1)*X2[2]-u*P(2,1)*X2[2])*w;
      A(ind1,7)=(P(0,1)*X2[3]-u*P(2,1)*X2[3])*w;
      A(ind1,8)= (P(0,2)*X2[0]-u*P(2,2)*X2[0])*w;
      A(ind1,9)= (P(0,2)*X2[1]-u*P(2,2)*X2[1])*w;
      A(ind1,10)=(P(0,2)*X2[2]-u*P(2,2)*X2[2])*w;
      A(ind1,11)=(P(0,2)*X2[3]-u*P(2,2)*X2[3])*w;
      if(!affine){
	A(ind1,12)=(P(0,3)*X2[0]-u*P(2,3)*X2[0])*w;
	A(ind1,13)=(P(0,3)*X2[1]-u*P(2,3)*X2[1])*w;
	A(ind1,14)=(P(0,3)*X2[2]-u*P(2,3)*X2[2])*w;
	A(ind1,15)=(P(0,3)*X2[3]-u*P(2,3)*X2[3])*w;
      }

      A(ind2,0)=(P(1,0)*X2[0]-v*P(2,0)*X2[0])*w;
      A(ind2,1)=(P(1,0)*X2[1]-v*P(2,0)*X2[1])*w;
      A(ind2,2)=(P(1,0)*X2[2]-v*P(2,0)*X2[2])*w;
      A(ind2,3)=(P(1,0)*X2[3]-v*P(2,0)*X2[3])*w;
      A(ind2,4)=(P(1,1)*X2[0]-v*P(2,1)*X2[0])*w;
      A(ind2,5)=(P(1,1)*X2[1]-v*P(2,1)*X2[1])*w;
      A(ind2,6)=(P(1,1)*X2[2]-v*P(2,1)*X2[2])*w;
      A(ind2,7)=(P(1,1)*X2[3]-v*P(2,1)*X2[3])*w;
      A(ind2,8)= (P(1,2)*X2[0]-v*P(2,2)*X2[0])*w;
      A(ind2,9)= (P(1,2)*X2[1]-v*P(2,2)*X2[1])*w;
      A(ind2,10)=(P(1,2)*X2[2]-v*P(2,2)*X2[2])*w;
      A(ind2,11)=(P(1,2)*X2[3]-v*P(2,2)*X2[3])*w;
      if(!affine){
	A(ind2,12)=(P(1,3)*X2[0]-v*P(2,3)*X2[0])*w;
	A(ind2,13)=(P(1,3)*X2[1]-v*P(2,3)*X2[1])*w;
	A(ind2,14)=(P(1,3)*X2[2]-v*P(2,3)*X2[2])*w;
	A(ind2,15)=(P(1,3)*X2[3]-v*P(2,3)*X2[3])*w;
      }
    }

    Matrix U,S,V;
    
    //svd(A) and svd(A'*A) give same vector, just its quicker to do latter
    //A.Lsvd(U,S,V);
    (A.transpose()*A).Lsvd(U,S,V);
    Matrix sol=V.getColumn(V.n-1);
        
    Matrix res=(A*sol);
    //printf("residual on it %d: %f\n",its,res.dot(res));
    
    //reshape H
    for(int i=0;i<(affine?3:4);i++){
      for(int j=0;j<4;j++){
	H(i,j)=sol[i*4+j];
      }
    }
  }
}
//assuming bottom row is [0,0,0,1]
Matrix closestEuclideanHomography(Matrix & H,
				  Matrix & R,Matrix & t,double & sc)
{
  //H.printMatlab("H");
  R=H.submatrix(0,0,3,3);
  Matrix U,S,V;
  R.Lsvd(U,S,V);
  R=U*V.transpose();
  sc=(S[0]+S[1]+S[2])/3;
  t=Matrix(3,1);
  t[0]=H(0,3)/sc;
  t[1]=H(1,3)/sc;
  t[2]=H(2,3)/sc;
  //printf("%f %f %f\n",S[0],S[1],S[2]);
  Matrix A=Matrix::eye(4,4);
  A.setSubmatrix(0,0,R);
  A.setSubmatrix(0,3,t);
  A=A*sc;
  A(3,3)=1;
  //A.printMatlab("A");
  return A;
}

vector<int> alignRobust(Matrix & H,vector<Matrix> & Ps,
			Matrix & X,Matrix & x,vector<int> & cams,
			int euclidean)
{
  printf("in align robust (euclidean:%d)\n",euclidean);
  const double thresh=4;
  vector<int> bestinliers;
  int bestset[8];
  int minpts=euclidean?6:8;

  double maxtime=0;
  for(int its=0;its<100;its++){
    if(its%20==0)printf("%d ",its);
    Matrix Hcur;
    int cset[8];
    unique_set(cset,X.n,minpts);

    Matrix Xset(4,minpts);
    Matrix xset(2,minpts);
    vector<int> camset;
    vector<int> curinliers;
    for(int j=0;j<minpts;j++){
      Xset(0,j)=X(0,cset[j]);
      Xset(1,j)=X(1,cset[j]);
      Xset(2,j)=X(2,cset[j]);
      Xset(3,j)=X(3,cset[j]);

      xset(0,j)=x(0,cset[j]);
      xset(1,j)=x(1,cset[j]);
      camset.push_back(cams[cset[j]]);
    }
    nacb::StopWatch swatch;
    swatch.start();
    alignImageBased(Hcur,Ps,Xset,xset,camset,euclidean); 
    double tm=(double)swatch;
  
    maxtime=std::max(tm,maxtime);

    if(euclidean){
      Matrix R,t;
      double sc;
      Hcur=closestEuclideanHomography(Hcur,R,t,sc);
      align_refine_euc(R,t,sc,Ps,Xset,xset,camset);
      Hcur=Matrix::eye(4,4);
      Hcur.setSubmatrix(0,0,R);
      Hcur.setSubmatrix(0,3,t);
      Hcur=Matrix::scale(sc,sc,sc)*Hcur;
    }
    else
      align_refine(Hcur,Ps,Xset,xset,camset);

    vector<Matrix> PHs;
    for(int i=0;i<(int)Ps.size();i++)
      PHs.push_back(Ps[i]*Hcur);

    for(int i=0;i<X.n;i++){
      Matrix & P=PHs[cams[i]];
      Matrix pp=projectPoints(P,X.getColumn(i));
      Vec2f diff(pp[0]-x(0,i),pp[1]-x(1,i));
    
      if(diff.len()<thresh){
	curinliers.push_back(i);
      }
    }
    if(curinliers.size()>bestinliers.size()){
      bestinliers=curinliers;
      memcpy(bestset,cset,sizeof(bestset));
      H=Hcur;
      //break if we found solution that agrees with all X
      if((int)bestinliers.size()==X.n)break;
    }
  }
  printf("\nmaxtime %f\n",maxtime);
  if((int)bestinliers.size()<minpts){
    printf("there are less than %d inliers\n",minpts);
    bestinliers.clear();
    return bestinliers;
  }
  maxtime=0;

  for(int its=0;its<2;its++){
    if((int)bestinliers.size()<minpts){
      printf("there are less than %d inliers\n",minpts);
      bestinliers.clear();
      return bestinliers;
    }

    Matrix Xbest(4,bestinliers.size());
    Matrix xbest(2,bestinliers.size());

    vector<int> camsbest;
    for(int i=0;i<(int)bestinliers.size();i++){
      Xbest(0,i)=X(0,bestinliers[i]);
      Xbest(1,i)=X(1,bestinliers[i]);
      Xbest(2,i)=X(2,bestinliers[i]);
      Xbest(3,i)=X(3,bestinliers[i]);
      
      xbest(0,i)=x(0,bestinliers[i]);
      xbest(1,i)=x(1,bestinliers[i]);
      camsbest.push_back(cams[bestinliers[i]]);
    }
    nacb::StopWatch swatch;
    
    Matrix Hnew;
    alignImageBased(Hnew,Ps,Xbest,xbest,camsbest);

    double tm=(double)swatch;
    maxtime=std::max(tm,maxtime);
    
    if(euclidean){
      Matrix R,t;
      double sc;
      Hnew=closestEuclideanHomography(Hnew,R,t,sc);
      align_refine_euc(R,t,sc,Ps,Xbest,xbest,camsbest);
      Hnew=Matrix::eye(4,4);
      Hnew.setSubmatrix(0,0,R);
      Hnew.setSubmatrix(0,3,t);
      Hnew=Matrix::scale(sc,sc,sc)*Hnew;
    }
    else
      align_refine(Hnew,Ps,Xbest,xbest,camsbest);

    //Get the inliers again with the better fit of the H matrix
    vector<Matrix> PHs;
    for(int i=0;i<(int)Ps.size();i++)
      PHs.push_back(Ps[i]*Hnew);
    
    vector<int> newbest;
    for(int i=0;i<X.n;i++){
      Matrix & P=PHs[cams[i]];
      Matrix pp=projectPoints(P,X.getColumn(i));
      Vec2f diff(pp[0]-x(0,i),pp[1]-x(1,i));
    
      if(diff.len()<thresh)newbest.push_back(i);
    }
    if(newbest.size()>bestinliers.size()){
      bestinliers=newbest;
      H=Hnew;
    }
    else break;
  }
  printf(" align_robust: inliers %ld/%d (final refine %f)\n",bestinliers.size(),X.n,maxtime);

  return bestinliers;
}

namespace TriOpt{
  int ncams;
  vector<Matrix> Ps;
  vector<Matrix> xs;
  int cnt=0;
  int jind;
  double g_err;
  void triangulate_err(const int & m,const int & n,
		       double * x,double * fvec,
		       int & info){
    Matrix X(4,1);
    X[0]=x[0];
    X[1]=x[1];
    X[2]=x[2];
    X[3]=1;

    g_err=0;
    for(int i=0;i<ncams;i++){
      Matrix pp=Ps[i]*X;
      fvec[2*i  ]=pp[0]/pp[2]-xs[i](0,jind);
      fvec[2*i+1]=pp[1]/pp[2]-xs[i](1,jind);	
      g_err+=fvec[2*i]*fvec[2*i];
      g_err+=fvec[2*i+1]*fvec[2*i+1];
    }
    if(isnan(g_err)||isinf(g_err))info=-1;
    g_err/=ncams;
  }
  void triangulate_an_err(const int & m,const int & n,
		       double * x,double * fvec,
		       double * fjac,const int & ldfjac,
		       int & info){
    Matrix X(4,1);
    X[0]=x[0];
    X[1]=x[1];
    X[2]=x[2];
    X[3]=1;
    if(info==1){
      g_err=0;
      for(int i=0;i<ncams;i++){
	Matrix pp=Ps[i]*X;
	fvec[2*i  ]=pp[0]/pp[2]-xs[i](0,jind);
	fvec[2*i+1]=pp[1]/pp[2]-xs[i](1,jind);	
	g_err+=fvec[2*i]*fvec[2*i];
	g_err+=fvec[2*i+1]*fvec[2*i+1];
      }
      g_err/=ncams;
      if(isnan(g_err)||isinf(g_err))info=-1;
    }
    else if(info==2){
      for(int i=0;i<ncams;i++){
	Matrix pp=Ps[i]*X;
	Matrix & aj=Ps[i];

	double imx=pp[0];
	double imy=pp[1];
	double imz=pp[2];

	double imzz=imz*imz;

	fjac[2*i         ]=aj[0]/imz-aj[8]*imx/imzz;
	fjac[2*i+  ldfjac]=aj[1]/imz-aj[9]*imx/imzz;
	fjac[2*i+2*ldfjac]=aj[2]/imz-aj[10]*imx/imzz;
	
	fjac[2*i+1        ]=aj[4]/imz-aj[8]*imy/imzz;
	fjac[2*i+1 +ldfjac]=aj[5]/imz-aj[9]*imy/imzz;
	fjac[2*i+1+2*ldfjac]=aj[6]/imz-aj[10]*imy/imzz;
      }
    }
  }
};

/** Perform the non-linear minimization of the triangulation
 *  of a set of points.
 * \param Ps a vector of camera matrices
 * \param xs a vector of 2 x npts matrices that store the observations
 * \param X  an input/output matrix of size 4 x npts, used as the 
 *           starting position
 * \return the average reprojection error of the points, and modifies
 *         the X matrix.
 */
double triangulate_opt(vector<Matrix> & Ps,
		       vector<Matrix> & xs,
		       Matrix & X)
{
  assert(xs[0].n==X.n);
  TriOpt::ncams=Ps.size();
  TriOpt::Ps=Ps;
  TriOpt::xs=xs;
  TriOpt::cnt=0;
  TriOpt::jind=0;
  Matrix x(3,1);
  Matrix fvec(2*xs.size(),1);
  double total_err=0;

  for(int j=0;j<X.n;j++){
    TriOpt::jind=j;

    x[0]=X(0,j)/X(3,j);
    x[1]=X(1,j)/X(3,j);
    x[2]=X(2,j)/X(3,j);
    
    int info;
    //TriOpt::triangulate_err(2*xs.size(),3,x.data,fvec.data,info);
    //printf("before tri-opt %f\n",sqrt(TriOpt::g_err));
    
    Matrix::lsqnonlin(TriOpt::triangulate_an_err,x,fvec);
    
    //get the error
    TriOpt::triangulate_err(2*xs.size(),3,x.data,fvec.data,info);
    
    X(0,j)=x[0];
    X(1,j)=x[1];
    X(2,j)=x[2];
    X(3,j)=1;
    
    total_err+=TriOpt::g_err;
  }
  return total_err/X.n;
}

vector<int>  triangulate_robust(vector<Matrix> & Ps,
				vector<Matrix> & xs,
				Matrix & X)
{
  const double thresh=2;
  const double thresh2=thresh*thresh;
  vector<int> bestinliers;

  for(int i=0;i<20;i++){
    int curset[2];
    unique_set(curset,xs.size(),2);
    vector<Matrix> Ps_set;
    vector<Matrix> xs_set;
    Ps_set.push_back(Ps[curset[0]]);
    Ps_set.push_back(Ps[curset[1]]);
    xs_set.push_back(xs[curset[0]]);
    xs_set.push_back(xs[curset[1]]);
    
    Matrix X=triangulate(Ps_set,xs_set);
    triangulate_opt(Ps_set,xs_set,X);

    vector<int> inliers;
    for(int i=0;i<(int)Ps.size();i++){
      Matrix d=projectPoints(Ps[i],X)-xs[i];
      if(d.dot(d)<thresh2){
	inliers.push_back(i);
      }
    }
    if(inliers.size()>bestinliers.size()){
      bestinliers=inliers;
      if(bestinliers.size()==xs.size())break;
    }
  }
  if(bestinliers.size()<2){
    bestinliers.clear();
    X=Matrix(4,1);
    X[0]=0;
    X[1]=0;
    X[2]=0;    
    X[3]=1;
    return bestinliers;
  }
  vector<Matrix> Ps_set;
  vector<Matrix> xs_set;
  for(int i=0;i<(int)bestinliers.size();i++){
    Ps_set.push_back(Ps[bestinliers[i]]);
    xs_set.push_back(xs[bestinliers[i]]);
  }
  X=triangulate(Ps_set,xs_set);
  triangulate_opt(Ps_set,xs_set,X);
  return bestinliers;
}

/** Compute a normalizing transform for a set of points.
 * \param p an 3xn vector of points (3rd row assumed to be 1)
 * \return a transformation that takes the points to the origin
 *         and scales them to be unit length.
 */
nacb::Matrix normalizingTransform(nacb::Matrix & p)
{
  nacb::Vec3d mean(0,0,0);
  double len=0;
  for(int i=0;i<p.n;i++)
    mean+=Vec3d(p(0,i),p(1,i),1);
  mean/=p.n;
  for(int i=0;i<p.n;i++){
    len+=(Vec3d(p(0,i),p(1,i),1)-mean).len();
  }
  len/=p.n;
  len+=1e-10;
  Matrix sc=Matrix::eye(3,3);
  Matrix tr=Matrix::eye(3,3);
  tr(0,2)=-mean.x;
  tr(1,2)=-mean.y;

  sc(0,0)=1.0/len;
  sc(1,1)=1.0/len;

  return sc*tr;
}


/** Create the skew symmetric matrix that corresponds to the 
 *  cross-product matrix of the vector v.
 * \param v the vector argument (assuming it is 3x1 or 1x3)
 * \return the 3x3 skew symmetric matrix
 */
Matrix skewSymmetric(const Matrix & v){
  Matrix E(3,3);
  E.setAll(0);
  E(0,1)=-v.get(2);
  E(0,2)= v.get(1);
  E(1,0)= v.get(2);
  E(1,2)=-v.get(0);
  E(2,0)=-v.get(1);
  E(2,1)= v.get(0);
  return E;
}

/** Convert projection matrices (first one must be eye(3,4)) into
 *  the corresponding trifocal tensor. From Hartley and Zisserman
 *  (pg. 366).
 * \param Ps the vector of projection matrices (Ps.size()==3)
 * \param T the output array containing the tensor (must be 3 matrices)
 */
void projection2trifocal(vector<Matrix> & Ps,Matrix T[])
{
  for(int k=0;k<3;k++){
    T[k]=Matrix(3,3);
    T[k].setAll(0);
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	T[i](j,k)=Ps[1](j,i)*Ps[2](k,3)-Ps[1](j,3)*Ps[2](k,i);
      }
    }
  }
}

/** Convert trifocal tensor into canonical projection matrices.
 *  From Hartley and Zisserman (pg.369).
 * \param T the array of matrices corresponding to the trifocal
 *        tensor
 * \return A vector of 3 projection matrices, the first in canonical
 *         form
 */
vector<Matrix> trifocal2projection(Matrix T[]){
  Matrix U,S,V,V1(3,3),V2(3,3);
  V2.setAll(0);
  V1.setAll(0);
  //Compute the right and left null spaces of T_i,
  //and place in cols of V2 and V1 respectively
  for(int k=0;k<3;k++){
    T[k].Lsvd(U,S,V);
    V2.setSubmatrix(k,0,V.getColumn(2).transpose());
    V1.setSubmatrix(k,0,U.getColumn(2).transpose());
  }
  //e1 is the right null space of V1 and
  //e2 is the right null space of V2
  V1.Lsvd(U,S,V);
  Matrix e1=V.getColumn(2);
  V2.Lsvd(U,S,V);
  Matrix e2=V.getColumn(2);
  Matrix e1x=skewSymmetric(e1);
  Matrix e2x=skewSymmetric(e2);
  
  Matrix F21(3,3);
  Matrix F31(3,3);
  for(int k=0;k<3;k++){
    Matrix v1=e1x*(T[k]*e2);
    Matrix v2=e2x*(T[k].transpose()*e1);
    F21.setSubmatrix(0,k,v1);
    F31.setSubmatrix(0,k,v2);
  }
  vector<Matrix> Ps;
  Ps.push_back(Matrix::eye(3,4));
  Ps.push_back(Matrix::eye(3,4));
  Ps.push_back(Matrix::eye(3,4));

  Ps[1].setSubmatrix(0,0,T[0]*e2);
  Ps[1].setSubmatrix(0,1,T[1]*e2);
  Ps[1].setSubmatrix(0,2,T[2]*e2);
  Ps[1].setSubmatrix(0,3,e1);
  
  Matrix M=e2*e2.transpose()-Matrix::eye(3,3);
  Ps[2].setSubmatrix(0,0,M*(T[0].transpose()*e1));
  Ps[2].setSubmatrix(0,1,M*(T[1].transpose()*e1));
  Ps[2].setSubmatrix(0,2,M*(T[2].transpose()*e1));
  Ps[2].setSubmatrix(0,3,e2);
  
  Ps[1]*=(1.0/sqrt(Ps[1].dot(Ps[1])));
  Ps[2]*=(1.0/sqrt(Ps[2].dot(Ps[2])));
  return Ps;
}

/** Compute the trifocal tensore using a linear algorithm for the
 *  correspondences in the 2 x npts matrices.
 * \param p1in the 2xnpts matrix of observations in 1st image
 * \param p2in the 2xnpts matrix of observations in 2st image
 * \param p3in the 2xnpts matrix of observations in 3st image
 * \return a vector of matrices corresponding to the canonical 
 *         projection matrices for the trifocal tensor.
 * \sa projection2trifocal trifocal2projection
 */
vector<Matrix> linearTrifocalTensor(Matrix & p1in,Matrix & p2in,Matrix & p3in)
{
  //First add an extra row of 1's to all the points
  int npts=p1in.n;
  Matrix p1(3,npts);
  Matrix p2(3,npts);
  Matrix p3(3,npts);
  for(int j=0;j<npts;j++){
    p1(0,j)=p1in(0,j);
    p1(1,j)=p1in(1,j);
    p1(2,j)=1;
    p2(0,j)=p2in(0,j);
    p2(1,j)=p2in(1,j);
    p2(2,j)=1;
    p3(0,j)=p3in(0,j);
    p3(1,j)=p3in(1,j);
    p3(2,j)=1;
  }
  //Compute and normalize the points
  Matrix H[3]={normalizingTransform(p1),
	       normalizingTransform(p2),
	       normalizingTransform(p3)};
  Matrix Hinv[3]={H[0].inverse(),
		  H[1].inverse(),
		  H[2].inverse()};
  p1=H[0]*p1;
  p2=H[1]*p2;
  p3=H[2]*p3;

  //The system of equations goes in A.
  //Tind gives the index into a linear vector of the
  //unknowns (basically assumes row-major)
  Matrix A(p1.n*4,27);
  int Tind[3][3][3];
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
	Tind[i][j][k]=i*9+j*3+k;
  
  A.setAll(0);
  //There are 4 constraints for each point correspondence
  int coni=0;
  for(int ind=0;ind<p1.n;ind++){
    for(int i=0;i<2;i++){
      for(int l=0;l<2;l++){
	for(int k=0;k<3;k++){
	  A(coni,Tind[k][2][2])= p1(k,ind)*p2(i,ind)*p3(l,ind);
	  A(coni,Tind[k][i][2])=-p1(k,ind)*p3(l,ind);
	  A(coni,Tind[k][2][l])=-p1(k,ind)*p2(i,ind);
	  A(coni,Tind[k][i][l])= p1(k,ind);
	}
	coni++;
      }
    }
  }
  //Get the solution using SVD and last column of V
  Matrix U,S,V;
  A.Lsvd(U,S,V);
  Matrix T[3]={Matrix(3,3),Matrix(3,3),Matrix(3,3)};
  Matrix t=V.getColumn(26);

  //Plug the solution into the T matrices
  for(int k=0;k<3;k++)
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	T[k](i,j)=t[Tind[k][i][j]];
  //printf("got T\n");

  //Now we must recover the unnormalized transformation.
  //First get the epipoles from the tensor.
  Matrix V2(3,3);
  Matrix V1(3,3);
  V2.setAll(0);
  V1.setAll(0);
  for(int k=0;k<3;k++){
    T[k].Lsvd(U,S,V);
    V2.setSubmatrix(k,0,V.getColumn(2).transpose());
    (T[k].transpose()).Lsvd(U,S,V);
    V1.setSubmatrix(k,0,V.getColumn(2).transpose());
  }
  V1.Lsvd(U,S,V);
  Matrix e1=V.getColumn(2);
  V2.Lsvd(U,S,V);
  Matrix e2=V.getColumn(2);
  
  //Now create the other system of equations that maps
  //the 18 unknowns of the 3x3 submatrices aij bij into
  //the 27 unknowns of the trifocal tensor
  Matrix E(27,18);
  E.setAll(0);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	E(Tind[i][j][k],  3*j+i)= e2[k];//a(j,i)
	E(Tind[i][j][k],9+3*k+i)=-e1[j];//b(k,i)
      }
    }
  }
  //minimize ||A*t||=||A*E*a||
  //solution is t=E*a
  E.Lsvd(U,S,V);
  Matrix Up=U.submatrix(0,0,27,18);
  
  Matrix AUp=A*Up;
  AUp.Lsvd(U,S,V);
  Matrix xp=V.getColumn(V.n-1);
  Matrix x=Up*xp;
  
  //x is the new trifocal tensor, put it back into T
  for(int k=0;k<3;k++)
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	T[k](i,j)=x[Tind[k][i][j]];
  
  //Now get the unnormalized solution
  Matrix Tnew[3]={Matrix(3,3),Matrix(3,3),Matrix(3,3)};
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	double res=0;
	for(int r=0;r<3;r++)
	  for(int s=0;s<3;s++)
	    for(int t=0;t<3;t++){
	      double v=H[0](r,i)*Hinv[1](j,s)*Hinv[2](k,t)*T[r](s,t);
	      res+=v;
	    }
	Tnew[i](j,k)=res;
      }
    }
  }
  //Convert the unnormalized tensor into projection matrices and return
  vector<Matrix> Ps=trifocal2projection(Tnew);
  projection2trifocal(Ps,Tnew);
  Ps=trifocal2projection(Tnew);
  return Ps;
}

vector<Matrix> getTrifocalTensorLinear(Matrix & p1in,Matrix & p2in,Matrix & p3in,
				       Matrix & p1eval,//all the points to evaluate the 3 soln's
				       Matrix & p2eval,
				       Matrix & p3eval,
				       int * bestscore_out)
{
  vector<Matrix> Ps=linearTrifocalTensor(p1in,p2in,p3in);
  
  vector<Matrix> pseval;
  
  pseval.push_back(p1eval);
  pseval.push_back(p2eval);
  pseval.push_back(p3eval);
  
  //Check this solution
  Matrix Xeval=triangulate(Ps,pseval);
  triangulate_opt(Ps,pseval,Xeval);
  
  int score=0;
  
  vector<int> inliers;
  Matrix dtot(Xeval.n,1);
  dtot.setAll(0);
  for(int j=0;j<3;j++){
    Matrix pp=projectPoints(Ps[j],Xeval);
    Matrix diff=pseval[j]-pp;
    for(int k=0;k<diff.n;k++){
      double d=diff(0,k)*diff(0,k)
	+diff(1,k)*diff(1,k);
      d=sqrt(d);
      dtot[k]+=d;
    }
  }
  dtot*=(1.0/3);
  for(int k=0;k<dtot.m;k++){
    if(dtot[k]<2){
      inliers.push_back(k);
      score++;
    }
  }
  if(inliers.size()>=12){
    Matrix p1i(2,inliers.size());
    Matrix p2i(2,inliers.size());
    Matrix p3i(2,inliers.size());
    for(int i=0;i<(int)inliers.size();i++){
      p1i(0,i)=p1eval(0,inliers[i]);
      p1i(1,i)=p1eval(1,inliers[i]);
      p2i(0,i)=p2eval(0,inliers[i]);
      p2i(1,i)=p2eval(1,inliers[i]);
      p3i(0,i)=p3eval(0,inliers[i]);
      p3i(1,i)=p3eval(1,inliers[i]);
    }
    Ps=linearTrifocalTensor(p1i,p2i,p3i);
    Xeval=triangulate(Ps,pseval);
    double avgerr=triangulate_opt(Ps,pseval,Xeval);
    
    if(avgerr<1000)
      bundleAdjustSparse(Ps,pseval,Xeval);
    score=0;
    dtot.setAll(0);
    for(int j=0;j<3;j++){
      Matrix pp=projectPoints(Ps[j],Xeval);
      Matrix diff=pseval[j]-pp;
      for(int k=0;k<diff.n;k++){
	double d=diff(0,k)*diff(0,k)
	  +diff(1,k)*diff(1,k);
	d=sqrt(d);
	dtot[k]+=d;
      }
    }
    dtot*=(1.0/3);
    inliers.clear();
    for(int k=0;k<dtot.m;k++){
      if(dtot[k]<2){
	inliers.push_back(k);
	score++;
      }
    }
  }
  if(bestscore_out)*bestscore_out=score;
  Ps.push_back(Xeval);
  return Ps;
}

/** Recover the trifocal tensor using a minimal parameterization.  This is from
 *  a paper by Torr.  The input points, p1in,p2in,p3in must have only 6 entries
 *  corresponding to 6 3D points in general form.  This function requires
 *  other observations in order to pick the best of 3 solutions to a cubic
 *  equation.  This are given as input in p1eval, p2eval and p3eval.
 *
 * \param p1in the 2x6 observation matrix of points in image 1
 * \param p2in the 2x6 observation matrix of points in image 2
 * \param p3in the 2x6 observation matrix of points in image 3
 * \param p1eval more point observations in image 1
 * \param p2eval more point observations in image 2
 * \param p3eval more point observations in image 3
 * \param bestscore_out is the number of inliers out of all the points in p*eval
 * \param returns the 3 projection matrices of the cameras
 */
vector<Matrix>  getTrilinearTensorMinimal(Matrix & p1in,Matrix & p2in,Matrix & p3in,
					  Matrix & p1eval,//all the points to evaluate the 3 soln's
					  Matrix & p2eval,
					  Matrix & p3eval,
					  int * bestscore_out)
{
  //add extra row of 1's and rectify the points (maps some observations to infinity)
  Matrix p1(3,6);
  Matrix p2(3,6);
  Matrix p3(3,6);
  for(int j=0;j<6;j++){
    p1(0,j)=p1in(0,j);
    p1(1,j)=p1in(1,j);
    p1(2,j)=1;
    p2(0,j)=p2in(0,j);
    p2(1,j)=p2in(1,j);
    p2(2,j)=1;
    p3(0,j)=p3in(0,j);
    p3(1,j)=p3in(1,j);
    p3(2,j)=1;
  }
  Matrix Bs[]={rectifyPoints(p1),
	       rectifyPoints(p2),
	       rectifyPoints(p3)};

  Matrix p[3]={Bs[0]*p1,
	       Bs[1]*p2,Bs[2]*p3};

  Matrix A(3,5);
  for(int i=0;i<3;i++){
    A(i,0)=-p[i](0,4)*p[i](1,5)+p[i](0,4)*p[i](2,5);
    A(i,1)= p[i](0,5)*p[i](1,4)-p[i](1,4)*p[i](2,5);
    A(i,2)=-p[i](0,5)*p[i](2,4)+p[i](1,5)*p[i](2,4);
    A(i,3)=-p[i](0,4)*p[i](2,5)+p[i](1,4)*p[i](2,5);
    A(i,4)= p[i](0,4)*p[i](1,5)-p[i](1,5)*p[i](2,4);
  }
  Matrix U,s,V;
  A.Lsvd(U,s,V);
  //t1 t2 t5 - t2 t3 t5 - t2 t4 t5 =
  // t1 t3 t4 - t2 t3 t4 -t3 t4 t5
  double a=0,b=0,c=0,d=0;
  Matrix t1(5,1),t2(5,1);
  for(int i=0;i<5;i++){
    t1[i]=V(i,3);
    t2[i]=V(i,4);
  }
  //t1.printMatlab("t1");
  //t2.printMatlab("t2");
  /*
    (a1+x b1)*(a2+x b2)*(a5+x b5)
    a1 a2 a5 + (a2 b1 +a1 b2) a5 x +b1 b2 a5 x^2
    a1 a2 x b5 + (a2 b1 +a1 b2) b5 x^2 +b1 b2 b5 x^3=
    a1 a2 a5+((a2 b1 +a1 b2) a5 + a1 a2 b5) x + (b1 b2 a5 + (a2 b1+a1 b2)b5)x^2
    +b1 b2 b5 x^3
  */
#define CUBE_TERM(sign,i,j,k) a+=sign*(t2[i]*t2[j]*t2[k]);\
  b+=sign*(t2[i]*t2[j]*t1[k]+(t1[j]*t2[i]+t1[i]*t2[j])*t2[k]);\
    c+=sign*((t1[j]*t2[i]+t1[i]*t2[j])*t1[k]+t1[i]*t1[j]*t2[k]);\
      d+=sign*(t1[i]*t1[j]*t1[k]);

    CUBE_TERM(+1.0,0,1,4);
    CUBE_TERM(-1.0,1,2,4);
    CUBE_TERM(-1.0,1,3,4);
    CUBE_TERM(-1.0,0,2,3);
    CUBE_TERM(+1.0,1,2,3);
    CUBE_TERM(+1.0,2,3,4);
    //printf("%f %f %f %f\n",a,b,c,d);
    double x[3];
    int nroots=cuberoots(a,b,c,d,x[0],x[1],x[2]);

    vector<Matrix> pseval;
    vector<Matrix> bestPs;
    Matrix bestX(4,1);
    int bestscore=0;

    pseval.push_back(p1eval);
    pseval.push_back(p2eval);
    pseval.push_back(p3eval);

    //printf("solution:\n");
    for(int i=0;i<nroots;i++){
      //printf("%f(%f) ",x[i],cubeeval(a,b,c,d,x[i]));
      Matrix t=t2*x[i]+t1;
      double diff=
	t[0]*t[1]*t[4]-t[1]*t[2]*t[4]-t[1]*t[3]*t[4]-
	(t[0]*t[2]*t[3]-t[1]*t[2]*t[3]-t[2]*t[3]*t[4]);
      //check the difference
      if(fabs(diff)>=1e-6)
	printf("NOTE: diff is %f!!!!!!! (should be 0, used to assert this)\n",diff);
      double X=(t[3]-t[4])/(t[1]-t[2]);
      double Y=t[3]/(t[0]-t[2]);
      double Z=t[4]/(t[0]-t[1]);
      double W=1;

      vector<Matrix> Ps;
      //if this happens the solution is probably bad,
      //so clamping will just make the residual high
      if(isinf(X)||isnan(X)||fabs(X)>10e13)X=10e13;
      if(isinf(Y)||isnan(Y)||fabs(Y)>10e13)Y=10e13;
      if(isinf(Z)||isnan(Z)||fabs(Z)>10e10)Z=10e13;
             
      for(int j=0;j<3;j++){
	Matrix A(4,4);
	A.setAll(0);
	A(0,0)= p[j](2,4);
	A(0,2)=-p[j](0,4);
	A(0,3)= p[j](2,4)-p[j](0,4);
	
	A(1,1)= p[j](2,4);
	A(1,2)=-p[j](1,4);
	A(1,3)= p[j](2,4)-p[j](1,4);
	
	A(2,0)= p[j](2,5)*X;
	A(2,2)=-p[j](0,5)*Z;
	A(2,3)= p[j](2,5)*W-p[j](0,5)*W;
	
	A(3,1)= p[j](2,5)*Y;
	A(3,2)=-p[j](1,5)*Z;
	A(3,3)= p[j](2,5)*W-p[j](1,5)*W;

	Matrix U,S,V;
	A.Lsvd(U,S,V);
	//S.printMatlab("S");
	Matrix cam=V.getColumn(3);

	 Matrix P=Matrix::eye(3,4);
	 P(0,0)=cam[0];
	 P(1,1)=cam[1];
	 P(2,2)=cam[2];
	 P(0,3)=cam[3];
	 P(1,3)=cam[3];
	 P(2,3)=cam[3];

	 P=Bs[j].inverse()*P;
	 //P.printMatlab("P");

	 Ps.push_back(P);


	 //check that the solution is correct
	 /*
	 Matrix X1(4,6);
	 X1.setAll(0);
	 X1(0,0)=1;
	 X1(1,1)=1;
	 X1(2,2)=1;
	 X1(3,3)=1;
	 X1(0,4)=1;
	 X1(1,4)=1;
	 X1(2,4)=1;
	 X1(3,4)=1;
	 X1(0,5)=X;
	 X1(1,5)=Y;
	 X1(2,5)=Z;
	 X1(3,5)=W;

	 Matrix pts=P*X1;
	 Matrix pts2(2,6);
	 for(int k=0;k<6;k++){
	   pts2(0,k)=pts(0,k)/pts(2,k);
	   pts2(1,k)=pts(1,k)/pts(2,k);
	 }
	 if(j==0)(pts2-p1in).printMatlab("p1diff");
	 if(j==1)(pts2-p2in).printMatlab("p2diff");
	 if(j==2)(pts2-p3in).printMatlab("p3diff");
	 */
       }
       Matrix T=Matrix::eye(4,4);
       Matrix P3x3inv=Ps[0].submatrix(0,0,3,3).inverse();
       T.setSubmatrix(0,0,P3x3inv);
       T.setSubmatrix(0,3,P3x3inv*Ps[0].submatrix(0,3,3,1)*-1);

       Matrix M=Matrix::eye(4,4);
       M(3,0)=-1;
       M(3,1)=-1;
       M(3,2)=-1;

       for(int j=0;j<3;j++)
	 Ps[j]=Ps[j]*T*M;

       //Trying to fix the bad P matrices, but this didn't work
       //Matrix Tfocal[3];
       //projection2trifocal(Ps,Tfocal);
       //Ps=trifocal2projection(Tfocal);
       
       //Check this solution
       Matrix Xeval=triangulate(Ps,pseval);
       //triangulate_opt(Ps,pseval,Xeval);

       Matrix dtot(Xeval.n,1);
       dtot.setAll(0);
       for(int j=0;j<3;j++){
	 Matrix pp=projectPoints(Ps[j],Xeval);
	 Matrix diff=pseval[j]-pp;
	 for(int k=0;k<diff.n;k++){
	   double d=diff(0,k)*diff(0,k)
	     +diff(1,k)*diff(1,k);
	   d=sqrt(d);
	   dtot[k]+=d;
	 }
	 //thisscore+=(diff.dot(diff));
       }
       dtot*=(1.0/3);
       int thisscore=0;
       for(int k=0;k<dtot.m;k++){
	 thisscore+=(dtot[k]<1);
       }
       if(thisscore>bestscore){
	 bestscore=thisscore;
	 bestPs=Ps;
	 bestX=Xeval;
       }
     }
     //printf("\n");
     //printf("bestscore is %f\n",bestscore/(pseval[0].n*3));
     if(bestscore_out){
       *bestscore_out=bestscore;
     }
     if(bestPs.size()==3 && bestPs[0].n==4){
       bestPs.push_back(bestX);
     }
     //bestPs[0].printMatlab("Ps[0]");
     return bestPs;
 }

/** Robustly fit a trifocal tensor to the image observations.
 * \param ps1 a 2xnpts matrix of observations in image 1
 * \param ps2 a 2xnpts matrix of observations in image 2
 * \param ps3 a 2xnpts matrix of observations in image 3
 * \param Xout the optional output is set to the inlying 3D points.
 * \param inds_out this optional output is set to be the indices
 *        into the observations of the inlying points
 *  \return The vector of projection matrices.
 */
vector<Matrix> robustTrilinear(Matrix & ps1,Matrix & ps2,Matrix & ps3,
			       Matrix * Xout,
			       vector<int> * inds_out)
{
  const int minpts=6;
  const double thresh=2;
  vector<Matrix> bestPs; 
  int bestscore=0;
  for(int its=0; its<1000; its++){
    if(its%200==0)printf("it:%d,bestscore:%d\n",its,bestscore);
    int curset[minpts];
    unique_set(curset,ps1.n,minpts);
    
    Matrix ps1use(2,minpts);
    Matrix ps2use(2,minpts);
    Matrix ps3use(2,minpts);
    
    for(int k=0;k<minpts;k++){
      ps1use(0,k)=ps1(0,curset[k]);
      ps1use(1,k)=ps1(1,curset[k]);
      ps2use(0,k)=ps2(0,curset[k]);
      ps2use(1,k)=ps2(1,curset[k]);
      ps3use(0,k)=ps3(0,curset[k]);
      ps3use(1,k)=ps3(1,curset[k]);
    }
    int score=0;
    vector<Matrix> Ps=getTrilinearTensorMinimal(ps1use,ps2use,ps3use,
    	ps1,ps2,ps3,&score);
    //vector<Matrix> Ps=getTrifocalTensorLinear(ps1use,ps2use,ps3use,
    //   ps1,ps2,ps3,&score);
    if(score>bestscore){
      bestscore=score;
      bestPs=Ps;
      //if we have found a solution with all inliers then break early
      if(bestscore==ps1.n)break;
      //printf("best score is %d\n",bestscore);
    }
  }
  const int print=0;
  vector<int> inds;//good indices
  double diff=0;
  Matrix X=bestPs[3];
  bestPs.pop_back();//remove points from the bestPs
  vector<Matrix> pseval;
  //X.printMatlab("X");
  
  pseval.push_back(ps1);
  pseval.push_back(ps2);
  pseval.push_back(ps3);
  
  Matrix dtot(X.n,1);
  Matrix zs(X.n,1);
  dtot.setAll(0);
  zs.setAll(0);

  for(int i=0;i<3;i++){
    int infront=0;

    Matrix PP=bestPs[i]*X;
    Matrix pp=projectPoints(bestPs[i],X);
    if(print)pp.printMatlab("pp");
    Matrix d=pseval[i]-pp;
    
    for(int k=0;k<d.n;k++){
      infront+=(PP(2,k)>0);
      dtot[k]+=sqrt(d(0,k)*d(0,k)+d(1,k)*d(1,k));
    }
    if(infront<d.n/2){
      bestPs[i]*=-1;
    }

    PP=bestPs[i]*X;
    for(int k=0;k<d.n;k++){
      zs[k]+=(PP(2,k)>0);
    }
  }
  dtot*=1.0/3;
  int cnt=0;
  for(int i=0;i<X.n;i++)
    if(dtot[i]<thresh && zs[i]==3)cnt++;
  Matrix Xgood(4,cnt);
  vector<Matrix> pobsgood;
  pobsgood.push_back(Matrix(2,cnt));
  pobsgood.push_back(Matrix(2,cnt));
  pobsgood.push_back(Matrix(2,cnt));
  
  cnt=0;
  for(int i=0;i<X.n;i++){
    if(dtot[i]<thresh && zs[i]==3){
      Xgood(0,cnt)=X(0,i);
      Xgood(1,cnt)=X(1,i);
      Xgood(2,cnt)=X(2,i);
      Xgood(3,cnt)=X(3,i);
      pobsgood[0](0,cnt)=ps1(0,i);
      pobsgood[0](1,cnt)=ps1(1,i);
      pobsgood[1](0,cnt)=ps2(0,i);
      pobsgood[1](1,cnt)=ps2(1,i);
      pobsgood[2](0,cnt)=ps3(0,i);
      pobsgood[2](1,cnt)=ps3(1,i);
      inds.push_back(i);
      cnt++;
    }
  }
  printf("%d/%d points are good\n",cnt,X.n);
  if(cnt<12){
    printf("need at least 12 points\n");
    return vector<Matrix>();
  }
  
  //Get new
  /*bestPs=linearTrifocalTensor(pobsgood[0],pobsgood[1],pobsgood[2]);
  Xgood=triangulate(bestPs,pobsgood);
  double avgres=triangulate_opt(bestPs,pobsgood,Xgood);
  if(isnan(avgres)||isinf(avgres)){
    printf("its nan somehow\n");
    return vector<Matrix>();
    }*/
  //end get and refine new
  /*
  //Compute T-form taking points to unit ball
  Matrix mn(3,1);
  mn.setAll(0);
  for(int j=0;j<Xgood.n;j++){
    mn[0]+=Xgood(0,j);
    mn[1]+=Xgood(1,j);
    mn[2]+=Xgood(2,j);
  }
  mn*=(1.0/Xgood.n);
  Matrix Xc(3,Xgood.n);
  for(int j=0;j<Xgood.n;j++){
    Xc(0,j)=Xgood(0,j);
    Xc(1,j)=Xgood(1,j);
    Xc(2,j)=Xgood(2,j);
  }
  Matrix covar=Xc*Xc.transpose();
  covar*=(1.0/Xgood.n);
  Matrix evec,eval,blah;
  covar.Lsvd(evec,eval,blah);
  Matrix sc=Matrix::eye(3,3);
  sc(0,0)=sqrt(eval[0])+1e-10;
  sc(1,1)=sqrt(eval[0])+1e-10;
  sc(2,2)=sqrt(eval[0])+1e-10;
  
  Matrix Hinv=Matrix::eye(4,4);
  Matrix T=Matrix::trans(-mn[0],-mn[1],-mn[2]);
  Hinv.setSubmatrix(0,0,(covar*sc).inverse());
  Hinv=Hinv*T;
  Matrix H=Hinv.inverse();
  bestPs[0]=bestPs[0]*H;
  bestPs[1]=bestPs[1]*H;
  bestPs[2]=bestPs[2]*H;
  Xgood=Hinv*Xgood;
  */

  //comparison of the sparse and full methods.
  if(0)
  {
    vector<Matrix> bestPs2;
    Matrix Xgood2=Xgood.copy();
    for(int i=0;i<(int)bestPs.size();i++)
      bestPs2.push_back(bestPs[i].copy());
    bundleAdjust(bestPs2,pobsgood,Xgood2);
  }
  bundleAdjustSparse(bestPs,pobsgood,Xgood);
  
  if(print)Xgood.printMatlab("X");
  printf("best score is %d, total diff:%f\n",bestscore,diff/(3*ps1.n));
  for(int i=0;i<3;i++){
    char name[256];
    sprintf(name,"P%d",i);
    if(print)bestPs[i].printMatlab(name);
  }
  
   //successful if half of the points are inliers
  if(inds.size()>ps1.n*0.5){
    if(Xout)*Xout=Xgood;
    if(inds_out)*inds_out=inds;
    return bestPs;
  }
  return vector<Matrix>();
}

Matrix getFundamentalMatrix(Matrix & p1in,Matrix & p2in){
  Matrix p1=p1in.copy();
  Matrix p2=p2in.copy();
  Vec2f x1c(0,0);
  Vec2f x2c(0,0);
  double x1s=0;
  double x2s=0;
  
  //Compute the centers
  for(int i=0;i<p1.n;i++){
    x1c.x+=p1(0,i);
    x1c.y+=p1(1,i);
  }
  x1c*=(1.0/((double)p1.n));

  for(int i=0;i<p2.n;i++){
    x2c.x+=p2(0,i);
    x2c.y+=p2(1,i);
  }
  x2c*=(1.0/((double)p2.n));

  //Compute the scales
  for(int i=0;i<p1.n;i++){
    p1(0,i)-=x1c.x;
    p1(1,i)-=x1c.y;
    
    double d=p1(0,i)*p1(0,i)+p1(1,i)*p1(1,i);
    x1s+=d;
  }
  x1s=sqrt(x1s/p1.n);

  for(int i=0;i<p2.n;i++){
    p2(0,i)-=x2c.x;
    p2(1,i)-=x2c.y;
    
    double d=p2(0,i)*p2(0,i)+p2(1,i)*p2(1,i);
    x2s+=d;
  }
  x2s=sqrt(x2s/p2.n);

  //scale down
  p1=p1*(1.0/x1s);
  p2=p2*(1.0/x2s);

  Matrix A(p1.n,9);
  for(int i=0;i<p1.n;i++){
    A(i,0)=p2(0,i)*p1(0,i);
    A(i,1)=p2(0,i)*p1(1,i);
    A(i,2)=p2(0,i);
    A(i,3)=p2(1,i)*p1(0,i);
    A(i,4)=p2(1,i)*p1(1,i);
    A(i,5)=p2(1,i);
    A(i,6)=p1(0,i);
    A(i,7)=p1(1,i);
    A(i,8)=1;
  }
  
  Matrix U,S,V;
  
  if(p1.n>1000){
    fprintf(stderr,"WARNING: defaulting to svd of A'*A...not gauranteed to work\n");
    fprintf(stderr,"WARNING(continued): modified sept 2008\n");
    (A.transpose()*A).Lsvd(U,S,V);
  }
  else
    A.Lsvd(U,S,V);
  
  Matrix F(3,3);
  for(int i=0;i<9;i++)
    F[i]=V(i,8);
  
  F.Lsvd(U,S,V);
  Matrix Smat=Matrix::eye(3,3);
  Smat(0,0)=S[0];
  Smat(1,1)=S[1];
  Smat(2,2)=0;//S[2];
  F=U*Smat*V.transpose();

  Matrix tr1=Matrix::eye(3,3);
  Matrix sc1=Matrix::eye(3,3);
  Matrix tr2=Matrix::eye(3,3);
  Matrix sc2=Matrix::eye(3,3);

  tr1(0,2)=-x1c.x;
  tr1(1,2)=-x1c.y;
  
  tr2(0,2)=-x2c.x;
  tr2(1,2)=-x2c.y;

  sc1(0,0)=1.0/x1s;
  sc1(1,1)=1.0/x1s;

  sc2(0,0)=1.0/x2s;
  sc2(1,1)=1.0/x2s;

  F=tr2.transpose()*sc2*F*sc1*tr1;

  /*
  for(int i=0;i<p1.n;i++){
    Matrix x1(3,1);
    Matrix x2(3,1);
    x1[0]=p1in(0,i);
    x1[1]=p1in(1,i);
    x1[2]=1;
    
    x2[0]=p2in(0,i);
    x2[1]=p2in(1,i);
    x2[2]=1;
    
    Matrix l1=F*x1;
    Matrix l2=x2.transpose()*F;
    
    l1*=(1.0/sqrt(l1[0]*l1[0]+l1[1]*l1[1]));
    l2*=(1.0/sqrt(l2[0]*l2[0]+l2[1]*l2[1]));
    
    Matrix d1=x2.transpose()*l1;
    Matrix d2=l2*x1;
    //printf("%f ",fabs(d1[0])+fabs(d2[0]));
  }
  //printf("\n");*/
  return F;
}

inline Vec3f fmat_norm_eline(Matrix & F,double x,double y){
  double lx=F.data[0]*x+F.data[1]*y+F.data[2];
  double ly=F.data[3]*x+F.data[4]*y+F.data[5];
  double lz=F.data[6]*x+F.data[7]*y+F.data[8];
  double len=(sqrt(lx*lx+ly*ly));
  return Vec3f(lx/len,ly/len,lz/len);
}

inline Vec3f fmat_norm_eline(double x,double y,Matrix & F){
  double lx=F.data[0]*x+F.data[3]*y+F.data[6];
  double ly=F.data[1]*x+F.data[4]*y+F.data[7];
  double lz=F.data[2]*x+F.data[5]*y+F.data[8];
  double len=(sqrt(lx*lx+ly*ly));
  return Vec3f(lx/len,ly/len,lz/len);
}
inline double dist_to_eline(Vec3f & eline,Vec2f & pt){
  return eline.x*pt.x+eline.y*pt.y+eline.z;
}

vector<int> fmatrix(vector<Feature *> & f1,	     
		    vector<Feature *> & f2, 
		    double ethresh,
		    double distanceLimit)
{
  const int bestk=3;       //was 2, with 1001 iterations
  const int nits = 3001; //was 2001
  bool useMedian = true;   //Use median score

  double worstscore=10e10;
  int * f1match=new int[f1.size()];
  int * f2match=new int[f2.size()];
  double * f1score=new double[f1.size()];
  double * f2score=new double[f2.size()];
  int * f1match2=new int[f1.size()];
  int * f2match2=new int[f2.size()];
  double * f1score2=new double[f1.size()];
  double * f2score2=new double[f2.size()];

  vector<IndexAndScore> * f1bestk=new vector<IndexAndScore>[f1.size()];
  vector<IndexAndScore> * f2bestk=new vector<IndexAndScore>[f2.size()];

  memset(f1match,0xFF,sizeof(int)*f1.size());
  memset(f2match,0xFF,sizeof(int)*f2.size());

  for(int i=0;i<(int)f1.size();i++)
    f1score[i]=worstscore;
  for(int i=0;i<(int)f2.size();i++)
    f2score[i]=worstscore;

  //get matching scores
  for(int i=0;i<(int)f1.size();i++){
    for(int j=0;j<(int)f2.size();j++){
      if((f1[i]->point - f2[j]->point).len()>=distanceLimit)
	continue;

      double sc=f1[i]->compare(f2[j]);

      IndexAndScore ias,jas;
      ias.index=j;
      ias.score=sc;
      jas.index=i;
      jas.score=sc;

      if(f1bestk[i].size()==0)
	f1bestk[i].push_back(ias);
      else if(f1bestk[i].size()<bestk ||
	      sc<f1bestk[i][f1bestk[i].size()-1].score){
	if(f1bestk[i].size()<bestk)f1bestk[i].push_back(ias);
	for(int k = 0; k<(int)f1bestk[i].size();k++){
	  if(sc<f1bestk[i][k].score){
	    for(int k2=(int)f1bestk[i].size()-1;k2>k;k2--)
	      f1bestk[i][k2]=f1bestk[i][k2-1];
	    f1bestk[i][k]=ias;
	    break;
	  }
	}
      }
      if(f2bestk[j].size()==0)
	f2bestk[j].push_back(jas);
      else if(f2bestk[j].size()<bestk ||
	      sc<f2bestk[j][f2bestk[j].size()-1].score){
	if(f2bestk[j].size()<bestk)f2bestk[j].push_back(jas);
	for(int k = 0; k<(int)f2bestk[j].size();k++){
	  if(sc<f2bestk[j][k].score){
	    for(int k2 = (int)f2bestk[j].size()-1;k2>k;k2--)
	      f2bestk[j][k2]=f2bestk[j][k2-1];
	    f2bestk[j][k]=jas;
	    break;
	  }
	}
      }
      
      if(sc<f1score[i]){
	f1score[i]=sc;
	f1match[i]=j;
      }
      if(sc<f2score[j]){
	f2score[j]=sc;
	f2match[j]=i;
      }
    }
  }

  //remove matches that are not forward-backward
  vector<int> good_matches;
  int good=0;
  int total=0;
  for(int i=0;i<(int)f1.size();i++){
    int match;
    if((match=f1match[i])>=0){
      total++;
      if(f2match[match]!=i){
	f1match[i]=-1;
      }
      else {
	good_matches.push_back(i);
	good++;
      }
    }
  }
  fprintf(stderr,"there are %d/%d\n",good,total);

  //can't fit if there is not enough good matches
  if(good<8){
    vector<int> res;
    return res;
  }

  //Now fit the fundamental matrix
  int bestset[8];
  double bestscore=10e10;

  vector<MatchAndScore> thescores;
  for(int its=0; its<nits; its++){
    if(its%1000==0)printf("its is %d (%f)\n",its,bestscore);
    int curset[8];

    for(int k=0; k<8; k++){
      int val=(int)(((double)rand())/((double)RAND_MAX)*(good_matches.size()));//%good_matches.size();
      int unique=1;
      for(int k2=0;k2<k;k2++){
	if(val==curset[k2]){
	  unique=0;
	  break;
	}
      }
      if(!unique){
	k--;
	continue;
      }
      curset[k]=val;
    }

    //fit f-matrix for curset
    Matrix p1(2,8);
    Matrix p2(2,8);
    for(int k=0;k<8;k++){
      int ind1=good_matches[curset[k]];
      int ind2=f1match[good_matches[curset[k]]];
      p1(0,k)=f1[ind1]->point.x;
      p1(1,k)=f1[ind1]->point.y;
      
      p2(0,k)=f2[ind2]->point.x;
      p2(1,k)=f2[ind2]->point.y;
    }
    Matrix F=getFundamentalMatrix(p1,p2);

    for(int i=0; i<(int)f1.size(); i++){
      double bestsc=10e10;
      int best=-1;      
      Vec3f l1 = fmat_norm_eline(F, f1[i]->point.x, f1[i]->point.y);

      //Used to maybe index out of bounds, now use f1bestk[i].size(), not bestk
      int kmax = f1bestk[i].size();

      for(int k=0; k<kmax; k++){
	int j=f1bestk[i][k].index;
	
	double d1=dist_to_eline(l1, f2[j]->point);
	if(fabs(d1)<bestsc){
	  best=j;
	  bestsc=fabs(d1);
	}
      }
      f1match2[i]=best;
      f1score2[i]=bestsc;
    }
    for(int j=0; j<(int)f2.size(); j++){
      double bestsc=worstscore;
      int best=-1;
      Vec3f l2 = fmat_norm_eline(f2[j]->point.x, f2[j]->point.y,F);
      
      int kmax = f2bestk[j].size();

      for(int k=0; k<kmax; k++){
	int i=f2bestk[j][k].index;
	
	double d2=dist_to_eline(l2,f1[i]->point);
	if(fabs(d2)<bestsc){
	  best=i;
	  bestsc=fabs(d2);
	}
      }
      f2match2[j]=best;
      f2score2[j]=bestsc;
    }
    vector<MatchAndScore> curscores;

    for(int i=0;i<(int)f1.size();i++){
      int j=f1match2[i];
      if(j>=0 && f2match2[j]==i){
	Vec3f l1 = fmat_norm_eline(F, f1[i]->point.x, f1[i]->point.y);
	Vec3f l2 = fmat_norm_eline(f2[j]->point.x, f2[j]->point.y, F);
	double d1 = dist_to_eline(l1, f2[j]->point);
	double d2 = dist_to_eline(l2, f1[i]->point);
	curscores.push_back(MatchAndScore(i,j,(fabs(d1)+fabs(d2))/2.0));
      }
    }
    std::sort(curscores.begin(),curscores.end());
    if(curscores.size()==0){
      printf("set is %d %d %d %d %d %d %d\n",
	     curset[0], curset[1], curset[2], curset[3],
             curset[4], curset[5], curset[6]);
	     
    }
    

    if(useMedian){
      if(curscores.size() && curscores[curscores.size()/2].score<bestscore){
	bestscore=curscores[curscores.size()/2].score;
	thescores=curscores;
	memcpy(bestset, curset, sizeof(curset));
      }
    }
    else {
      int inliers=0;
      
      for(int i=0; i<(int)f1.size(); i++){
	int j=f1match2[i];

	if(j>=0 && f2match2[j]==i){
	  Vec3f l1 = fmat_norm_eline(F, f1[i]->point.x, f1[i]->point.y);
	  Vec3f l2 = fmat_norm_eline(f2[j]->point.x, f2[j]->point.y, F);
	  double d1 = dist_to_eline(l1, f2[j]->point);
	  double d2 = dist_to_eline(l2, f1[i]->point);

	  curscores.push_back(MatchAndScore(i, j, (fabs(d1)+fabs(d2))/2.0));

	  if(fabs(d1)<ethresh && fabs(d2)<ethresh){
	    inliers--;
	  }
	}
      }
      //inliers is negative to be consistent with taking min of median score
      if(inliers<bestscore){
	bestscore = inliers;
	thescores = curscores;
	memcpy(bestset,curset,sizeof(curset));
      }
    }
  }

  vector<int> theinliers;
  vector<int> theinliers_match;
  for(int i=0;i<(int)thescores.size();i++){
    if(thescores[i].score<ethresh){
      theinliers.push_back(thescores[i].index1);
      theinliers_match.push_back(thescores[i].index2);
    }
  }
  printf("bestinliers: %ld/%ld (useMedian: %d)\n", theinliers.size(), good_matches.size(), useMedian);
  printf("match sizes: %ld,%ld\n", theinliers.size(), theinliers_match.size());

  Matrix m1(2, theinliers.size());
  Matrix m2(2, theinliers.size());
  
  for(int i=0;i<(int)theinliers.size();i++){
    m1(0,i)=f1[theinliers[i]]->point.x;
    m1(1,i)=f1[theinliers[i]]->point.y;

    m2(0,i)=f2[theinliers_match[i]]->point.x;
    m2(1,i)=f2[theinliers_match[i]]->point.y;
  }

  /*m1.printMatlab("m1");
  m2.printMatlab("m2");
  printf("hold on;\nfor i=1:size(m1,2),\nplot([m1(1,i),m2(1,i)],[m1(2,i),m2(2,i)]);\nend\n");
  printf("waitforbuttonpress;\n");
  */

  delete [] f1match;
  delete [] f2match;
  delete [] f1score;
  delete [] f2score;
  delete [] f1match2;
  delete [] f2match2;
  delete [] f1score2;
  delete [] f2score2;
  delete [] f1bestk;
  delete [] f2bestk;

  vector<int> results;
  
  //was 0.4, but makes more sense to be 1/2
  if(theinliers.size()>good_matches.size()*0.4){
    results = vector<int>(f1.size(), -1);

    for(int i=0;i<(int)theinliers.size();i++){
      results[theinliers[i]]=theinliers_match[i];
    }
  }
  return results;
}


/**getHomography XP=H*X
 *
 * normalize input points, and then use linear algorithm
 * to find the solution.
 * copied from calibrate.cc in svn/thesis/calibration
 */
Matrix getHomography(Matrix & X_,Matrix & XP_,int dlt)
{
  Matrix A(2*X_.n,8);
  Matrix b(2*X_.n,1);

  if( X_.m!=2 && X_.n>=4){
    static Matrix junk(1,1);
    printf("wrong dimensions for getHomography %s, %d\n",__FILE__,__LINE__);
    return junk;
  }
  Matrix X=X_.copy();
  Matrix XP=XP_.copy();
  Vec2f m(0,0);
  Vec2f mp(0,0);
  double s=0;
  double sp=0;
  //Compute mean 
  for(int j=0;j<X.n;j++){
    m.x+=X(0,j);
    m.y+=X(1,j);
    mp.x+=XP(0,j);
    mp.y+=XP(1,j);
  }
  m*=(1.0/(double)X.n);
  mp*=(1.0/(double)XP.n);
 
  //Subtract mean and compute scale factor
  for(int j=0;j<X.n;j++){
    X(0,j)-=m.x;
    X(1,j)-=m.y;
    XP(0,j)-=mp.x;
    XP(1,j)-=mp.y;
    s+=(sqrt(X(0,j)*X(0,j)+X(1,j)*X(1,j)));
    sp+=(sqrt(XP(0,j)*XP(0,j)+XP(1,j)*XP(1,j)));
  }
  double sqrt2=sqrt(2.0);
  s*=(sqrt2/(double)X.n);
  sp*=(sqrt2/(double)XP.n);
  s=1.0/s;
  sp=1.0/sp;
  X*=s;
  XP*=sp;

  Matrix H(3,3);
  if(dlt){
    A=Matrix(2*X_.n,9);
    for(int i=0;i<X.n;i++){
      double wp=1;
      double w =1;
      
      A(2*i,0)=0;
      A(2*i,1)=0;
      A(2*i,2)=0;
      A(2*i,3)=-wp*X(0,i);
      A(2*i,4)=-wp*X(1,i);
      A(2*i,5)=-wp*w;
      A(2*i,6)=XP(1,i)*X(0,i);
      A(2*i,7)=XP(1,i)*X(1,i);
      A(2*i,8)=XP(1,i)*w;
      A(2*i+1,0)=wp*X(0,i);
      A(2*i+1,1)=wp*X(1,i);
      A(2*i+1,2)=wp*w;
      A(2*i+1,3)=0;
      A(2*i+1,4)=0;
      A(2*i+1,5)=0;
      A(2*i+1,6)=-XP(0,i)*X(0,i);
      A(2*i+1,7)=-XP(0,i)*X(1,i);
      A(2*i+1,8)=-XP(0,i)*w;
    }
    Matrix u,s,v;
    A.Lsvd(u,s,v);
    H=v.getColumn(8);
    H.m=3;
    H.n=3;
  }
  else{
    //x'*z'=(a1*x+a2*y+a3) --> (a1*x0+a2*y0+a3)-x0'*z0'=0
    //y'*z'=(a4*x+a5*y+a6) --> (a4*x0+a5*y0+a6)-y0'*z0'=0
    //z'=(a7*x+a8*y+a9) 
    //(a1*x0+a2*y0+a3)-x0'*(a7*x+a8*y+a9)
    // = a1*x0+a2*y0+a3-x0'*a7*x-x0'*a8*y-x0' force a9==1
    for(int i=0;i<X.n;i++){
      A(2*i,0)=X(0,i);
      A(2*i,1)=X(1,i);
      A(2*i,2)=1;
      A(2*i,3)=0;
      A(2*i,4)=0;
      A(2*i,5)=0;
      A(2*i,6)=-XP(0,i)*X(0,i);
      A(2*i,7)=-XP(0,i)*X(1,i);
      b(2*i,0)=XP(0,i);
      A(2*i+1,0)=0;
      A(2*i+1,1)=0;
      A(2*i+1,2)=0;
      A(2*i+1,3)=X(0,i);
      A(2*i+1,4)=X(1,i);
      A(2*i+1,5)=1;
      A(2*i+1,6)=-XP(1,i)*X(0,i);
      A(2*i+1,7)=-XP(1,i)*X(1,i);
      b(2*i+1,0)=XP(1,i);
    }
    
    //Solve the equations
    Matrix a=Matrix::LlinLeastSq(A,b);
  
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	if(!(i==2 && j==2)){
	  H(i,j)=a[3*i+j];
	}
	else H(i,j)=1;
      }
    }
  }

  //Compute the associated scale and translation matrices
  Matrix Sp=Matrix::eye(3,3);
  Matrix Tp=Matrix::eye(3,3);
  Matrix S=Matrix::eye(3,3);
  Matrix T=Matrix::eye(3,3);
  Sp(0,0)=sp;
  Sp(1,1)=sp;
  Tp(0,2)=-mp.x;
  Tp(1,2)=-mp.y;
  S(0,0)=s;
  S(1,1)=s;
  T(0,2)=-m.x;
  T(1,2)=-m.y;
  //put the normalization in H
  //(Sp Tp)^-1 H S T 
  H=(Sp*Tp).inverse()*H*S*T;
  X=X_;
  XP=XP_;

  //Check
  Matrix X1s(3,X.n);
  X1s.setAll(1);
  for(int j=0;j<X.n;j++){
    X1s(0,j)=X(0,j);
    X1s(1,j)=X(1,j);
  }
  X1s=H*X1s;
  for(int j=0;j<X.n;j++){
    X1s(0,j)/=X1s(2,j);
    X1s(1,j)/=X1s(2,j);
  }
  Matrix diff(2,X.n);
  for(int j=0;j<X.n;j++){
    diff(0,j)=X1s(0,j)-XP(0,j);
    diff(1,j)=X1s(1,j)-XP(1,j);
  }
  //X.printMatlab("X");
  //XP.printMatlab("XP");
  //diff.printMatlab("diff");
  //H.printMatlab("H");
  //printf("HOMOGRAPHY ERROR: %f\n",
  // sqrt(diff.dot(diff)));
  return H;
}


double homographyModelError(Matrix & x1, Matrix & x2){
  double err=0;
  Matrix H=getHomography(x1,x2);
  Matrix Hinv=H.inverse();
  for(int k=0;k<x1.n;k++){
    Matrix x1_(3,1);
    Matrix x2_(3,1);
    x1_[0]=x1(0,k);
    x1_[1]=x1(1,k);
    x1_[2]=1;
    
    x2_[0]=x2(0,k);
    x2_[1]=x2(1,k);
    x2_[2]=1;
    Matrix x1p=H*x1_;
    x1p*=(1.0/x1p[2]);
    
    Matrix x2p=Hinv*x2_;
    x2p*=(1.0/x2p[2]);
    
    Matrix d1=x1p-x2_;
    Matrix d2=x2p-x1_;
    err+=(sqrt(d1.dot(d1))+sqrt(d2.dot(d2)))/2.0;
  }
  return err/x1.n;
}

Matrix projectiveFactorization(vector<Matrix> & obs,
			       vector<Matrix> & Ps){
  assert(obs.size());
  
  //make sure that all the observation vectors have
  //the same number of points
  int npts=obs[0].n;
  int ncams=obs.size();

  for(int i=0;i<(int)obs.size();i++)
    if(obs[i].n!=npts)
      return Matrix(1,1);
  
  Matrix A(3*ncams,npts);
  Matrix lams(ncams,npts);
  lams.setAll(1);
  
  vector<Matrix> Tinvs;
  vector<Matrix> uvs;
  for(int i=0;i<(int)obs.size();i++){
    Vec2f mn(0,0);
    double sc=0;
    for(int j=0;j<obs[i].n;j++){
      mn+=Vec2f(obs[i](0,j),obs[i](1,j));
    }
    mn*=(1.0/obs[i].n);
    Matrix T=Matrix::eye(3,3);
    
    Matrix ob(2,obs[i].n);
    for(int j=0;j<obs[i].n;j++){
      ob(0,j)=obs[i](0,j)-mn.x;
      ob(1,j)=obs[i](1,j)-mn.y;
      Vec2f n(ob(0,j),ob(1,j));
      sc+=n.len();
    }
    sc*=(1.0/obs[i].n);
    
    ob*=(1.0/sc);
   
    T(0,2)=mn.x;
    T(1,2)=mn.y;
    T(0,0)=sc;
    T(1,1)=sc;
    
    Tinvs.push_back(T);
    uvs.push_back(ob);
  }

  Matrix X,Psmat;
  const int maxits=20;
  for(int its=0;its<maxits;its++){
    //populate the A matrix
    for(int i=0;i<(int)obs.size();i++){
      for(int j=0;j<npts;j++){
	A(3*i  ,j)=uvs[i](0,j)*lams(i,j);
	A(3*i+1,j)=uvs[i](1,j)*lams(i,j);
	A(3*i+2,j)=lams(i,j);
      }
    }
    for(int re=0;re<3;re++){
      //rescale the columns and rows
      for(int j=0;j<npts;j++){
	double len=0;
	for(int i=0;i<A.m;i++)
	  len+=A(i,j)*A(i,j);
	len=sqrt(len);
	len=1.0/len;
	for(int i=0;i<A.m;i++)
	  A(i,j)*=len;
      }
      for(int i=0;i<(int)obs.size();i++){
	double len=0;
	for(int j=0;j<npts;j++){
	  len+=A(3*i,j)*A(3*i,j);
	  len+=A(3*i+1,j)*A(3*i+1,j);
	  len+=A(3*i+2,j)*A(3*i+2,j);
	}
	len=1.0/sqrt(len);
	
	for(int j=0;j<npts;j++){
	  A(3*i+0,j)*=len;
	  A(3*i+1,j)*=len;
	  A(3*i+2,j)*=len;
	  
	}
      }
    }
    Matrix U,S,V;
    A.Lsvd(U,S,V);
    Psmat=U.submatrix(0,0,U.m,4);
    X=V.submatrix(0,0,V.m,4);
    Matrix sc=Matrix::eye(4,4);
    for(int k=0;k<4;k++)
      sc(k,k)=sqrt(S[k]);

    Psmat=Psmat*sc;
    X=sc*X.transpose();
    
    for(int i=0;i<(int)obs.size();i++){
      for(int j=0;j<npts;j++){
	lams(i,j)=Psmat(3*i+2,0)*X(0,j)+Psmat(3*i+2,1)*X(1,j)+
	  Psmat(3*i+2,2)*X(2,j)+Psmat(3*i+2,3)*X(3,j);
      }
    }
  }
  Ps.clear();
  for(int i=0;i<(int)obs.size();i++){
    Ps.push_back(Tinvs[i]*Psmat.submatrix(3*i,0,3,4));
  }
  for(int i=0;i<X.n;i++){
    X(0,i)/=X(3,i);
    X(1,i)/=X(3,i);
    X(2,i)/=X(3,i);
    X(3,i)=1;
  }
  return X;
}


//Instantiate functions
template Vec3<float> triangulate(const Matrix & P1, const Matrix & P2,
				 const Vec2<float> & p1, const Vec2<float> & p2);
template Vec3<double> triangulate(const Matrix & P1, const Matrix & P2,
				  const Vec2<double> & p1, const Vec2<double> & p2);

#ifdef TEST_TRIFOCAL_ALG

double myrand(double min,double max)
{
  double r=((double)rand())/((double)RAND_MAX);
  return r*(max-min)+min;
}

int main(int ac,char * av[])
{
  Matrix K=Matrix::eye(3,3);
  K(0,0)=K(1,1)=500;
  K(0,2)=320;
  K(1,2)=240;
  
  vector<Matrix> Ps;
  vector<Matrix> xs;
  Matrix X=Matrix::random(4,20);
  for(int i=0;i<X.n;i++)
    X(3,i)=1;
  X=Matrix::trans(-0.5,-0.5,-0.5)*X;
  X.printMatlab("X");
  for(int k=0;k<3;k++){
    Ps.push_back(K*Matrix::eye(3,4)*
		 Matrix::trans(myrand(-1,1),myrand(-1,1),2)*
		 Matrix::rotx(myrand(-1,1))*
		 Matrix::roty(myrand(0,2*M_PI)));
    Ps[k]*=(1.0/sqrt(Ps[k].dot(Ps[k])));
  }
  Matrix Tr=Matrix::eye(4,4);
  Matrix P3x3inv=Ps[0].submatrix(0,0,3,3).inverse();
  Tr.setSubmatrix(0,0,P3x3inv);
  Tr.setSubmatrix(0,3,P3x3inv*Ps[0].submatrix(0,3,3,1)*-1);
  for(int i=0;i<(int)Ps.size();i++)Ps[i]=Ps[i]*Tr;
  for(int k=0;k<3;k++){
    xs.push_back(projectPoints(Ps[k],X));
    xs[k].printMatlab("x");
    for(int j=0;j<xs[k].n;j++){
      xs[k](0,j)+=myrand(-2,2);
      xs[k](1,j)+=myrand(-2,2);
    }
  }
  
  vector<Matrix> Ps3=robustTrilinear(xs[0],xs[1],xs[2]);
  vector<Matrix> Ps2=linearTrifocalTensor(xs[0],xs[1],xs[2]);
  /*Matrix T[3];
  projection2trifocal(Ps,T);
  vector<Matrix> Ps2=trifocal2projection(T);*/
  
  Matrix X3=triangulate(Ps3,xs);
  triangulate_opt(Ps3,xs,X3);
  
  Matrix X2=triangulate(Ps2,xs);
  double res=triangulate_opt(Ps2,xs,X2);
  printf("residual is %f\n",res);
  
  X2.printMatlab("X2");
  for(int i=0;i<3;i++){
    Matrix pp2=projectPoints(Ps2[i],X2);
    Matrix diff2=pp2-xs[i];
    Matrix pp3=projectPoints(Ps3[i],X3);
    Matrix diff3=pp3-xs[i];
    
    printf("%d %f %f\n",i,sqrt(diff2.dot(diff2)/X.n),
	   sqrt(diff3.dot(diff3)/X.n));
    
  }
}
#endif
