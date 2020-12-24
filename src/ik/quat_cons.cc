#include "quat_cons.h"
#include <nmath/matrix.h>
#include <nmath/vec3.h>
#include <lbfgs.h>
#include <iostream>

using namespace std;
using namespace nacb;

inline Quaternion quatXYZ(double * x){
  return Quaternion::rod(Vec3d(0,0,x[2])).times
    (Quaternion::rod(Vec3d(0,x[1],0)).times(Quaternion::rod(Vec3d(x[0],0,0))));
}
inline Quaternion quatXYZ_dx(double * x){
  return Quaternion::rod(Vec3d(0,0,x[2])).times
    (Quaternion::rod(Vec3d(0,x[1],0)).times
     (Quaternion::rod_dx(Vec3d(x[0],0,0))));
}
inline Quaternion quatXYZ_dy(double * x){
  return Quaternion::rod(Vec3d(0,0,x[2])).times
    (Quaternion::rod_dy(Vec3d(0,x[1],0)).times
     (Quaternion::rod(Vec3d(x[0],0,0))));
}
inline Quaternion quatXYZ_dz(double * x){
  return Quaternion::rod_dz(Vec3d(0,0,x[2])).times
    (Quaternion::rod(Vec3d(0,x[1],0)).times
     (Quaternion::rod(Vec3d(x[0],0,0))));
}

double quatDiffGradXYZ(int n,double * x,double * g,const void * qptr){
  Quaternion & q=*(Quaternion*)qptr;

  Quaternion c=quatXYZ(x);
  double da=c.a-q.a;
  Vec3d  dv=c.v-q.v;

  Quaternion qx=quatXYZ_dx(x);
  Quaternion qy=quatXYZ_dy(x);
  Quaternion qz=quatXYZ_dz(x);
    
  g[0]=2*da*qx.a+2*(dv.dot(qx.v));
  g[1]=2*da*qy.a+2*(dv.dot(qy.v));
  g[2]=2*da*qz.a+2*(dv.dot(qz.v));
  
  return da*da+dv.dot(dv);
}

double boundedClosestXYZ(const Quaternion & q,//should be const
		       double x[],
		       double l[],
		       double u[]){
  int nbd[3]={2,2,2};

  //clamp as lbfgs doesn't like stuff out of bounds
  for(int k=0;k<3;k++)x[k]=std::min(std::max(l[k],x[k]),u[k]);

  double res=lbfgs(quatDiffGradXYZ,3,x,l,u,nbd,&q);

  double xmid[3];
  for(int k=0;k<3;k++)xmid[k]=(l[k]+u[k])/2.0;
  double rmid=lbfgs(quatDiffGradXYZ,3,xmid,l,u,nbd,&q,1e4,1e-7);
  if(rmid<res){
    //printf("middle wins %f %f %g\n",rmid,res,res-rmid);
    memcpy(x,xmid,sizeof(double)*3);
    res=rmid;
  }
  return res;
}

#ifdef QUAT_CONS_MAIN
double rand_ang(double min,double max){
  return (double(rand())/double(RAND_MAX))*(max-min)+min;
}
int main(){
  Vec3d x=Vec3d(0,0,0)*(M_PI/180.0);
  Vec3d lg=Vec3d(-90,-90,-90)*(M_PI/180.0);
  Vec3d ug=Vec3d( 90 ,-80, -80)*(M_PI/180.0);

  Vec3d l=Vec3d(-60,-90,-90)*(M_PI/180.0);
  Vec3d u=Vec3d( 90 ,90, 90)*(M_PI/180.0);

  for(int i=0;i<10000;i++){
    //double rx=rand_ang(l.x,u.x);//
    double rx=rand_ang(lg.x,ug.x);
    double ry=rand_ang(lg.y,ug.y);
    double rz=rand_ang(lg.z,ug.z);
    char ind[3]={0,1,2};
    Vec3d d(rx,ry,rz);
    Quaternion des=quatXYZ(d.data);
    des.toEuler(ind,x.data);
    //x=x*0;
    boundedClosestXYZ(des,x.data,l.data,u.data);
    
    Quaternion res=quatXYZ(x.data);
    res.normalize();
    des.normalize();
    double dd=std::min(1.0,res.dot(des));
    double ac=acos(dd)*2*180/M_PI;
    if(ac>1e-3 && 0){
      cout << "res" << res << endl;
      cout << "des" << des << endl;

      printf("acos(dots)*2=%f  %f\n",acos(dd)*2*180/M_PI,
	     (fabs(x.data[0]-rx)+fabs(x.data[1]-ry)+
	     fabs(x.data[2]-rz))*180/M_PI); 
    }
  }
}
#endif
