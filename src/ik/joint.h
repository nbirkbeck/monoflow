#ifndef JOINT_H
#define JOINT_H

#pragma interface

#include <math.h>
#include <assert.h>
#include <GL/gl.h>

#include <nmath/matrix.h>
#include <nmath/vec3.h>
#include <nmath/mat4x4.h>
#include <nmath/quaternion.h>
#include <stdexcept>
#include <vector>
#include "quat_cons.h"

using nacb::Quaternion;

class Joint{
 public:
  bool useloc;      //!<Boolean for if we have translation dof (should be true for root only)
  bool lockloc;     //!<True if loction is set but fixed/locked (e.g., not a dof)  
  nacb::Vec3d      loc;   //!The animated location of the bone
  
  Joint(bool useloc=false,bool lockloc=false):loc(0.0,0.0,0.0){
    this->useloc=useloc;
    this->lockloc=lockloc;
  }
  virtual ~Joint(){ }
  virtual nacb::Vec3d getPoseLoc() const {
    return loc;
  }
  virtual Quaternion getPoseQuat() const {
    return Quaternion(1,nacb::Vec3d(0,0,0));
  }
  virtual void setPoseLoc(const nacb::Vec3d & loc){
    //FIXME: only allow if not locked?
    this->loc=loc;
  }
  virtual void setPoseQuat(const Quaternion & q){}
  virtual int applyDeriv(nacb::Mat4x4 & tform,int index,bool differential){
    if(useloc){
      if(lockloc)tform.trans_right(loc.x,loc.y,loc.z);
      else{
	if(index>=0 && index<3){
	  nacb::Vec4d col(tform(0,index),
			  tform(1,index),
			  tform(2,index),
			  tform(3,index));
	  tform.zero();
	  tform(0,3)=col.x;
	  tform(1,3)=col.y;
	  tform(2,3)=col.z;
	  tform(3,3)=col.w;
	}
	else tform.trans_right(loc.x,loc.y,loc.z);

	index-=3;//Only if !lockloc (so that nparms is same)
      }
    }
    return index;
  }
  virtual int getNParams() const {return (useloc&&!lockloc)?3:0;}
  virtual int getParams(double * p,bool differential) const {
    if(useloc && !lockloc){
      p[0]=loc.x;
      p[1]=loc.y;
      p[2]=loc.z;
    }
    return (useloc&&!lockloc)?3:0;
  }
  virtual int setParams(const double * p,bool differential){
    if(useloc && !lockloc){
      loc.x=p[0];
      loc.y=p[1];
      loc.z=p[2];
    }
    return (useloc&&!lockloc)?3:0;
  }

  virtual void glTransform() const {glTranslatef(loc.x,loc.y,loc.z);}
  virtual void apply(nacb::Mat4x4 & mat) const {
    mat.trans_right(loc.x,loc.y,loc.z);
  }
  virtual void writeBounds(FILE * file)const {}


  //Move this
  inline double clampAngle(double ang){
    while(ang<-4.0*M_PI)ang+=2.0*M_PI;
    while(ang> 4.0*M_PI)ang-=2.0*M_PI;
    return ang;
  }
};



class EulerJointXYZ: public Joint{
 public:
  double minx,maxx;//!< Joint limits for x-rotation
  double miny,maxy;//!< Joint limits for y-rotation
  double minz,maxz;//!< Joint limits for z-rotation
  
  double rx,ry,rz; //!<Rotation about x,y and z
  bool usex,usey,usez;

  EulerJointXYZ(bool usex=true,
		bool usey=true,
		bool usez=true,
		double minx=-M_PI*1000,double maxx=M_PI*1000,
		double miny=-M_PI*1000,double maxy=M_PI*1000,
		double minz=-M_PI*1000,double maxz=M_PI*1000):Joint(){
    this->usex=usex;
    this->usey=usey;
    this->usez=usez;

    this->minx=minx;
    this->miny=miny;
    this->minz=minz;

    this->maxx=maxx;
    this->maxy=maxy;
    this->maxz=maxz;

    rx=ry=rz=0;

    printf("ranges: %f %f   %f %f   %f %f\n",minx,maxx,miny,maxy,minz,maxz);
  }
  virtual Quaternion getPoseQuat() const{
    return 
      Quaternion::rod(nacb::Vec3d(0,0,usez?rz:0.0))*
      Quaternion::rod(nacb::Vec3d(0,usey?ry:0.0,0))*
      Quaternion::rod(nacb::Vec3d(usex?rx:0.0,0,0));
  }
  virtual void setPoseQuat(const Quaternion & q){
    if((int(usex)+int(usey)+int(usez))==1){
      rx=ry=rz=0;
      // There was a bug here (was using acos)
      if(usex)rx = atan2(q.v.x, q.a)*2.0;
      if(usey)ry = atan2(q.v.y, q.a)*2.0;
      if(usez)rz = atan2(q.v.z, q.a)*2.0;
      //printf("%d,%f %d,%f %d,%f  %f,%f,%f,%f\n", usex, rx, usey, ry, usez, rz,
      //     q.v.x, q.v.y, q.v.z, q.a);
      return;
    }
    char order[3];
    int k=0;
    double angles[3] = {0,0,0};
    double old[3]={rx,ry,rz};
    if(usex)order[k++]=0;
    if(usey)order[k++]=1;
    if(usez)order[k++]=2;
    
    if(!usex)order[k++]=0;
    if(!usey)order[k++]=1;
    if(!usez)order[k++]=2;

    q.toEuler(order,angles);
    
    //printf("order=[%d,%d,%d], angles=[%f,%f,%f]\n",order[0],order[1],order[2],
    //angles[0],angles[1],angles[2]);
    
    //unpack
    k=0;
    if(usex)rx=angles[k++];
    else rx=0;

    if(usey)ry=angles[k++];
    else ry=0;

    if(usez)rz=angles[k++];
    else rz=0;

    if(isnan(angles[0]) || isinf(angles[0]) ||
       isnan(angles[1]) || isinf(angles[1]) ||
       isnan(angles[2]) || isinf(angles[2])){
      Quaternion qc = q;
      fprintf(stderr, "WARNING: nan angles: %f %f %f nangles=%d %f>1 = %d\n", angles[0], angles[1], angles[2], k, qc.len(), qc.len()>1);
      if(isnan(rx) || isinf(rx))rx = old[0];
      if(isnan(ry) || isinf(ry))ry = old[1];
      if(isnan(rz) || isinf(rz))rz = old[2];
    }
       
    //printf("%f %f %f\n",rx,ry,rz);
  }
  virtual int applyDeriv(nacb::Mat4x4 & tform,int index,bool differential){
    index=Joint::applyDeriv(tform,index,differential);

    double rz_val=rz;
    double ry_val=ry;
    double rx_val=rx;
    
    if(differential && usex && usey && usez){
      if(usez)tform.rotz_right(rz);
      if(usey)tform.roty_right(ry);
      if(usex)tform.rotx_right(rx);
      
      rz_val=ry_val=rx_val=0;
    }

    if(usez){
      if(index==0)tform.dz_right(rz_val);
      else tform.rotz_right(rz_val);
      index--;
    }
    if(usey){
      if(index==0)tform.dy_right(ry_val);
      else tform.roty_right(ry_val);
      index--;
    }
    if(usex){
      if(index==0)tform.dx_right(rx_val);
      else tform.rotx_right(rx_val);
      index--;
    }
    return index;
  }
  virtual int getNParams() const {
    return Joint::getNParams()+
      int(usex)+int(usey)+int(usez);}
  virtual int getParams(double * p,bool differential) const {
    int n=Joint::getParams(p,differential);
    if(differential && usex && usey && usez){
      if(usez)p[n++]=0;
      if(usey)p[n++]=0;
      if(usex)p[n++]=0;
    }
    else{
      if(usez)p[n++]=rz;
      if(usey)p[n++]=ry;
      if(usex)p[n++]=rx;
    }
    return n;
  }
  virtual int setParams(const double * p,bool differential){
    int n=Joint::setParams(p,differential);
    if(differential && usex && usey && usez){
      double x[3]={rx,ry,rz};
      double l[3]={minx,miny,minz};
      double u[3]={maxx,maxy,maxz};

      nacb::Mat4x4 R;
      R.to_eye();
      R.rotz_right(rz);
      R.roty_right(ry);
      R.rotx_right(rx);


      double dz=p[n++];
      double dy=p[n++];
      double dx=p[n++];

      nacb::Mat4x4 Dz=nacb::Mat4x4::dz(0.0);
      nacb::Mat4x4 Dy=nacb::Mat4x4::dy(0.0);
      nacb::Mat4x4 Dx=nacb::Mat4x4::dx(0.0);
      nacb::Mat4x4 X=(Dz*dz)+(Dy*dy)+(Dx*dx);

      nacb::Mat4x4 X2=X*X;
      nacb::Mat4x4 X3=X2*X;
      nacb::Mat4x4 D=nacb::Mat4x4::eye()+X+(X2)*0.5+X3*(1.0/6.0);

      //cerr << "params: "<< p[n-3] << "," << p[n-2] << "," << p[n-1] << endl;
      R.mult_right(D);

      Quaternion q=Quaternion::fromMatrix(R);
      
      /*
      Quaternion q=
	Quaternion::rod(nacb::Vec3d(0,0,rz))*
	Quaternion::rod(nacb::Vec3d(0,ry,0))*
	Quaternion::rod(nacb::Vec3d(rx,0,0))*
	Quaternion::rod(nacb::Vec3d(0,0,dz))*
	Quaternion::rod(nacb::Vec3d(0,dy,0))*
	Quaternion::rod(nacb::Vec3d(dx,0,0)); 
      */
      
      boundedClosestXYZ(q,x,l,u);

      nacb::Vec3d d(x[0]-rx,x[1]-ry,x[2]-rz);
      if(d.len()>0.1){
	fprintf(stderr,"Norm changed quite a bit (%f)\n",d.len());
	fprintf(stderr,"Before %f %f %f   after %f %f %f\n",
		rx,ry,rz,x[0],x[1],x[2]);
      }
      
      ///char ind[3]={0,1,2};
      ///q.toEuler(ind,x);
      rz=clampAngle(x[2]);
      ry=clampAngle(x[1]);
      rx=clampAngle(x[0]);
      /*
      rz=std::min(std::max(rz,minz),maxz);
      ry=std::min(std::max(ry,miny),maxy);
      rx=std::min(std::max(rx,minx),maxx);
      */
    }
    else{
      if(usez){
	rz=clampAngle(p[n++]);
	rz=std::min(std::max(rz,minz),maxz);
      }
      if(usey){
	ry=clampAngle(p[n++]);
	ry=std::min(std::max(ry,miny),maxy);
      }
      if(usex){
	rx=clampAngle(p[n++]);
	rx=std::min(std::max(rx,minx),maxx);
      }
    }
    return n;
  }

  virtual void glTransform() const {
    Joint::glTransform();
    if(usez)glRotatef(this->rz*180.0/M_PI,0,0,1);
    if(usey)glRotatef(this->ry*180.0/M_PI,0,1,0);
    if(usex)glRotatef(this->rx*180.0/M_PI,1,0,0);
  }

  virtual void apply(nacb::Mat4x4 & mat) const {
    Joint::apply(mat);
    if(usez)mat.rotz_right(rz);
    if(usey)mat.roty_right(ry);
    if(usex)mat.rotx_right(rx);
  }

  virtual void writeBounds(FILE * file)const {
    Joint::writeBounds(file);
    double rad2deg=180.0/M_PI;
    if(usex)fprintf(file,"rx(%g,%g) ",minx*rad2deg,maxx*rad2deg);
    if(usey)fprintf(file,"ry(%g,%g) ",miny*rad2deg,maxy*rad2deg);
    if(usez)fprintf(file,"rz(%g,%g) ",minz*rad2deg,maxz*rad2deg);
  }
};






class RodJoint : public Joint {
 public:
  nacb::Vec3d      rod;   //!The angle/axis representation of rotation
  bool lockrot;

  RodJoint(const nacb::Vec3d & r = nacb::Vec3d(0.,0.,0.), 
	   bool useloc = false,
	   bool lockr = false): Joint(useloc), rod(r), lockrot(lockr) { }
  virtual ~RodJoint(){ }

  virtual Quaternion getPoseQuat() const {
    return Quaternion::rod(rod);
  }

  virtual void setPoseQuat(const Quaternion & q){
    double len = acos(q.a)*2; //a = cos(theta/2) => theta = acos(a)*2

    rod = q.v;
    rod.normalize();

    rod *= len;
  }

  // This function could benefit from caching.
  virtual int applyDeriv(nacb::Mat4x4 & tform, int index, bool differential){
    index = Joint::applyDeriv(tform, index, differential);

    int nrot = (lockrot ? 0 : 3);
    if (index < 0 || index >= nrot) return index - nrot;

    nacb::Mat4x4    skew;
    nacb::Mat4x4   dskew;
    
    nacb::Vec3d axisUn = rod;
    nacb::Vec3d axis = axisUn;
    nacb::Vec3d dangle(axis.x, axis.y, axis.z);
        
    double angle = std::max(axis.normalize(), 1e-5);
    double sa = sin(angle), ca = cos(angle);

    dangle *= (1.0/angle);

    memset(skew.data, 0, sizeof(nacb::Mat4x4));
    memset(dskew.data, 0, sizeof(nacb::Mat4x4));

    // Axis is normalized.
    skew(1, 2) = -axis.x;
    skew(2, 1) =  axis.x;

    skew(2, 0) = -axis.y;
    skew(0, 2) =  axis.y;

    skew(0, 1) = -axis.z;
    skew(1, 0) =  axis.z;
       
    /*
      (d/dx)rx/angle ==> (angle * ddx rx - rx * ddx angle)/(angle*angle)
    */
    
    dskew(1, 2) = -(angle * (index==0) - axisUn.x * dangle.data[index])/(angle*angle);
    dskew(2, 1) =  (angle * (index==0) - axisUn.x * dangle.data[index])/(angle*angle);
    
    dskew(2, 0) = -(angle * (index==1) - axisUn.y * dangle.data[index])/(angle*angle);
    dskew(0, 2) =  (angle * (index==1) - axisUn.y * dangle.data[index])/(angle*angle);
    
    dskew(0, 1) = -(angle * (index==2) - axisUn.z * dangle.data[index])/(angle*angle);
    dskew(1, 0) =  (angle * (index==2) - axisUn.z * dangle.data[index])/(angle*angle);
 
    //nacb::Mat4x4 rot = nacb::Mat4x4::eye() + skew * sa + (skew * skew)*(1 - ca);

    //memset(derivs, 0, sizeof(Mat4x4)*6);

    tform *= (dskew*sa + skew*(ca*dangle[index]) 
	      + (skew*dskew + dskew*skew)*(1-ca) + skew*skew*(sa*dangle[index]));
    return index-3;
  }

  virtual int getNParams() const {
    int ret = Joint::getNParams() + (lockrot ? 0 : 3);
    fprintf(stderr, "Getting rod nparms: %d\n", ret);
    return ret;
  }

  virtual int getParams(double * p,bool differential) const {
    int pi = Joint::getParams(p, differential);
    
    if (!lockrot) {
      p[pi    ] = rod.x;
      p[pi + 1] = rod.y;
      p[pi + 2] = rod.z;
    }
    return pi + (lockrot ? 0 : 3);
  }

  virtual int setParams(const double * p,bool differential){
    int pi = Joint::setParams(p, differential);

    if (!lockrot) {
      rod.x = p[pi + 0];
      rod.y = p[pi + 1];
      rod.z = p[pi + 2];
    }

    return pi + (lockrot ? 0 : 3);
  }

  virtual void glTransform() const {
    Joint::glTransform();
    Quaternion::rod(rod).glRotate();
  }

  virtual void apply(nacb::Mat4x4 & mat) const {
    Joint::apply(mat);
    mat *= nacb::Mat4x4(nacb::Vec3d(0,0,0), 
			nacb::Quaternion::rod(rod));
  }

  virtual void writeBounds(FILE * file) const {
    //FIXME: is this okay?
    fprintf(file, "aa() ");
  }
};




/** \brief This a more general bone type that has a number of 
           sequential degrees of freedom.  Some things are not 
	   completely implemented.
	    
  This class was created for the MIT joints, so it has a rotation
  around -90 and 90 before the animation.  I couldn't figure out the
  mapping to get these to work without, and it doesn't matter anyway.
*/
class EulerJoint: public Joint{
public:

  /** \brief Simple class representing a degree of freedom.
   */
  class Dof {
  public:
    Dof(int _coord = 0, double _value = 0, double _minValue = -1000, double _maxValue = 1000):
      coord(_coord), value(_value), minValue(_minValue), maxValue(_maxValue) { }
    
    void clamp() {
      value = std::max(minValue, std::min(maxValue, value));
    }

    nacb::Quaternion getQuaternion() const {
      double rx = (coord == 0) ? value : 0.0;
      double ry = (coord == 1) ? value : 0.0;
      double rz = (coord == 2) ? value : 0.0;
      return Quaternion::rod(nacb::Vec3d(rx, ry, rz));
    }

    nacb::Mat4x4 getMatrix() const {
      if (coord == 0) return nacb::Mat4x4::rotx(value);
      if (coord == 1) return nacb::Mat4x4::roty(value);
      if (coord == 2) return nacb::Mat4x4::rotz(value);
      return nacb::Mat4x4::eye();
    }

    void applyRight(nacb::Mat4x4& tform) const {
      if (coord == 0) tform.rotx_right(value);
      if (coord == 1) tform.roty_right(value);
      if (coord == 2) tform.rotz_right(value);
    }

    void applyDerivRight(nacb::Mat4x4& tform) const {
      if (coord == 0) tform.dx_right(value);
      if (coord == 1) tform.dy_right(value);
      if (coord == 2) tform.dz_right(value);
    }

    int coord;      //!< The index (x = 0, y = 1, z = 2)
    double value;   //!< Current value of the dof.
    double minValue;//!< The minValue
    double maxValue;//!< The max value.
  };

  std::vector<Dof> dof;

  EulerJoint(const std::vector<Dof>& _dof):Joint(), dof(_dof) {}

  //!\brief Get the pose quat (it does not include the extra rotation around x)
  virtual Quaternion getPoseQuat() const{
    nacb::Quaternion q;
    //Quaternion::rod(nacb::Vec3d(-M_PI/2.0, 0, 0));
    for (int i=0; i<(int)dof.size(); i++)
      q = q*dof[i].getQuaternion();
    //return q * Quaternion::rod(nacb::Vec3d(M_PI/2.0, 0, 0));
    return q;
  }

  //!\brief Set the pose quat (uses new code for getting the 2DOF cases)
  virtual void setPoseQuat(const Quaternion & q){
    if (!dof.size()) return;

    if(dof.size() == 1) {
      if (dof[0].coord == 0)
	dof[0].value = atan2(q.v.x, q.a)*2.0;

      if (dof[0].coord == 1)
	dof[0].value = atan2(q.v.y, q.a)*2.0;

      if (dof[0].coord == 2)
	dof[0].value = atan2(q.v.z, q.a)*2.0;
      return;
    }

    // Treat these guys separately; this helps!
    if (dof.size() == 2) {
      nacb::Mat3x3 m;
      q.getMatrix(m);

      bool debug = false;

      if (debug) {
	std::cout << "------Pose quat----" << std::endl;
	printf("%d %d\n", dof[0].coord, dof[1].coord);
	std::cout << getPoseQuat() << std::endl;
	std::cout << q << std::endl;
	printf("%f %f\n", dof[0].value, dof[1].value);
      }

      // Rx * Rz
      if (dof[0].coord == 0 && dof[1].coord == 2) {
	dof[0].value = atan2(-m(1, 2), m(2, 2));
	dof[1].value = atan2(-m(0, 1), m(0, 0));
      }
      
      else if (dof[0].coord == 0 && dof[1].coord == 1) {
	dof[0].value = atan2(m(2, 1), m(1, 1));
	dof[1].value = atan2(m(0, 2), m(0, 0));
      }

      else if (dof[0].coord == 1 && dof[1].coord == 2) {
	dof[0].value = atan2(m(0, 2), m(2, 2));
	dof[1].value = atan2(m(1, 0), m(1, 1));
      }
      else {
	// Use matlab syms to determine the above angles...
	printf("Unknown pair:%d %d\n", dof[0].coord, dof[1].coord);
      }

      if (debug) {
	printf("%d %d\n", dof[0].coord, dof[1].coord);
	std::cout << getPoseQuat() << std::endl;
      }

      return;
    }

    char order[3];
    int k=0;
    double angles[3] = {0,0,0};
    double old[3]={0,0,0};
    bool used[3] = {0, 0, 0};

    for (int i=int(dof.size() - 1); i >= 0; i--, k++) {
      order[k] = dof[i].coord;
      used[dof[i].coord] = 1;
      old[i] = dof[i].value;
    }

    // Order has to have all 3 angles (in the case of 2)
    for (int i=0; i<3; i++) {
      if (!used[i]) {
	order[k++] = i;
      }
    }

    bool debug = false;// dof.size() == 2;

    if (debug) {
      std::cout << "------Pose quat----" << std::endl;
      printf("%d %d %d\n", order[0], order[1], order[2]);
      std::cout << getPoseQuat() << std::endl;
      std::cout << q << std::endl;
      printf("old:  %f %f %f\n", old[0], old[1], old[2]);
    }

    q.toEuler(order,angles);

    if (dof.size() == 2)
      std::cout << angles[2] << std::endl;

    //unpack
    k = 0;
    for (int i=int(dof.size() - 1); i >= 0; i--, k++) {
      dof[i].value = angles[k];
    }
    if (debug) {
      printf("new:  %f %f %f\n", angles[2], angles[1], angles[0]);
      std::cout << getPoseQuat() << std::endl;
    }
    if(isnan(angles[0]) || isinf(angles[0]) ||
       isnan(angles[1]) || isinf(angles[1]) ||
       isnan(angles[2]) || isinf(angles[2])){
      Quaternion qc = q;
      fprintf(stderr, "WARNING: nan angles: %f %f %f nangles=%d %f>1 = %d\n", 
	      angles[0], angles[1], angles[2], k, qc.len(), qc.len()>1);

      for (int i=0; i<(int)dof.size(); i++) {
	if (isnan(dof[i].value))
	  dof[i].value = old[i];
      }
    }       
    //printf("%f %f %f\n",rx,ry,rz);
  }

  virtual int applyDeriv(nacb::Mat4x4 & tform,int index,bool differential){
    index=Joint::applyDeriv(tform,index,differential);
    
    // Make sure to apply the other rotations around x
    tform.rotx_right(-M_PI/2.0);
    for (int i=0; i<(int)dof.size(); i++) {
      if (index == 0)
	dof[i].applyDerivRight(tform);
      else
	dof[i].applyRight(tform);
      index--;
    }
    tform.rotx_right( M_PI/2.0);
    return index;
  }

  virtual int getNParams() const {
    return Joint::getNParams()+dof.size();
  }

  virtual int getParams(double * p,bool differential) const {
    int n=Joint::getParams(p,differential);
    for (int i=0; i<(int)dof.size(); i++, n++) {
      p[n] = dof[i].value;
    }
    return n;
  }

  virtual int setParams(const double * p,bool differential){
    int n=Joint::setParams(p,differential);
    for (int i=0; i<(int)dof.size(); i++, n++) {
      dof[i].value = p[n];
      dof[i].clamp();
    }
    return n;
  }

  virtual void glTransform() const {
    Joint::glTransform();
    glRotatef(-90, 1, 0, 0);
    for (int i=0; i<(int)dof.size(); i++) {
      dof[i].getQuaternion().glRotate();
    }
    glRotatef( 90, 1, 0, 0);
  }

  virtual void apply(nacb::Mat4x4 & mat) const {
    Joint::apply(mat);
    mat.rotx_right(-M_PI/2.0);    

    for (int i=0; i<(int)dof.size(); i++) 
      dof[i].applyRight(mat);

    mat.rotx_right( M_PI/2.0);
  }

  virtual void writeBounds(FILE * file)const {
    Joint::writeBounds(file);
    double rad2deg=180.0/M_PI;

    // Notice that the direction here is flipped..
    fprintf(file, "mit ");
    for (int i=0; i<(int)dof.size(); i++) {
      fprintf(file, "r%c(%g,%g) vx(%g) ", dof[i].coord + 'x', 
	      dof[i].minValue * rad2deg, dof[i].maxValue * rad2deg, dof[i].value * rad2deg); 
    }
  }
};


#endif
