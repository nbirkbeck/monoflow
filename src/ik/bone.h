#pragma interface

#ifndef BONE_H
#define BONE_H

#include <dgemm.h>
#include <assert.h>

#include <vector>
#include <list>
#include <GL/gl.h>
#include <GL/glu.h>
#include <map>

#include <nmath/mat4x4.h>
#include <nimage/image.h>
#include <nmath/matrix.h>
#include <nmath/vec3.h>
#include <nmath/vec4.h>
#include <nmath/quaternion.h>
#include <nmath/vec2.h>
#include <nmath/mat3x3.h>

#include <nmisc/timer.h>
#include "globals.h"
#include "joint.h"

using nacb::Matrix;
using nacb::Vec3d;
using nacb::Vec3f;
using nacb::Vec4d;
using nacb::Vec4f;
using nacb::Mat3x3;
using nacb::Mat4x4;
using nacb::Quaternion;

using std::vector;
using std::pair;
using std::string;
using nacb::Timer;
using nacb::StopWatch;

template <class obj,typename arg> 
struct my_binder{
public:
  obj & o;
  void (obj::*ptr)(const arg &); 
  my_binder(obj & o_in,
	    void (obj::*ptr_in)(const arg &)):o(o_in),ptr(ptr_in){
  }
  void operator()(const arg & a) const{
    ((&o)->*(ptr))(a);
  };
};

template <typename O,typename t>
my_binder<O,t> my_bind(O & o,void (O::*ptr)(const t &)){
  return my_binder<O,t>(o,ptr);
}

class Bone{
public:
  enum select_t{
    NotSelected=0x0,
    HeadSelected=0x1,
    TipSelected=0x2,
    BoneSelected=0x3
  }selected;

  //the indices into Constraints
  enum {
    TipLocCons,
    HeadLocCons,
    RotXCons,
    RotYCons,
    RotZCons,
    NumCons
  };

  static double cons_scales[NumCons];

  class Constraints{
  public:
    Vec4d data[NumCons];
    bool has_cons[NumCons];
    int n;
    Constraints(){
      n=0;
      memset(has_cons,0,sizeof(has_cons));
    }

    void write(FILE * f){
      unsigned short mask = 0;
      
      for(int i=0; i<NumCons; i++)
	mask |= (int((has_cons[i])>0) << i);
      
      fwrite(&mask, sizeof(short), 1, f);

      for(int i=0; i<NumCons; i++)
	if(has_cons[i])fwrite(data[i].data, sizeof(Vec4d), 1, f);
    }

    void read(FILE * f){
      unsigned short mask = 0;
      
      n = 0;
      
      size_t r = fread(&mask, sizeof(short), 1, f);
      if (r != 1) {
        fprintf(stderr, "%s:%d error reading mask\n", __FILE__, __LINE__);
      }
      for(int i=0; i<NumCons; i++){
	has_cons[i] = (mask>>i)&0x1;
	if(has_cons[i]){
	  r = fread(data[i].data, sizeof(Vec4d), 1, f);
          if (!r) {
            fprintf(stderr, "%s:%d error reading mask\n", __FILE__, __LINE__);
          }
	  n++;
	}
      }
    }

    void set(int k, const Vec4d & c){
      n += !(has_cons[k]);
      has_cons[k] = true;
      data[k] = c;
    }

    const Vec4d & get(int index) const{
      return data[index];
    }
    Vec4d & operator[](int index){
      return data[index];
    }
    void trans(double x,double y,double z){
      for(int k=0; k<NumCons; k++){
	if(has_cons[k]){
	  if(data[k].w==1){
	    data[k].x+=x;
	    data[k].y+=y;
	    data[k].z+=z;
	  }
	}
      }
    }
    void dx(double t){
      for(int k=0; k<NumCons; k++)if(has_cons[k])data[k].dx(t);
    }
    void rotx(double t){
      for(int k=0; k<NumCons; k++)if(has_cons[k])data[k].rotx(t);
    }
    void dy(double t){
      for(int k=0; k<NumCons;k++)if(has_cons[k])data[k].dy(t);
    }
    void roty(double t){
      for(int k=0; k<NumCons;k++)if(has_cons[k])data[k].roty(t);
    }
    void dz(double t){
      for(int k=0; k<NumCons;k++)if(has_cons[k])data[k].dz(t);
    }
    void rotz(double t){
      for(int k=0; k<NumCons;k++)if(has_cons[k])data[k].rotz(t);
    }
    void mult_left(const Mat4x4 & mat){
      for(int k=0; k<NumCons;k++)if(has_cons[k])data[k]=mat*data[k];
    }
    bool operator==(const Constraints & c) const{
      return memcmp(this,&c,sizeof(Constraints));
    }
  };
  Constraints cons;

  typedef pair<Bone*,uint32_t>       bone_affects_t;
  typedef Bone *                       constraint_t;
  typedef vector<constraint_t>     constraint_set_t;
  typedef pair<constraint_t,vector<Bone*> > chain_t;

  static bool drawBounds;
  static int updated;
public:
  double stiff;    //!<The stiffness (not used!)
  double len;      //!<The length of the bone
  vector<Bone *> kids;//!< List of children of this bone
  Bone * parent;   //!< Bones parent
  string name;     //!< Name of the bone
  int      id;     //!< Unique name

  int  cons_depth;//!< The depth that the constraint goes up the chain (0 means all the way)

  //Matrix lmat;//Local matrix
  //Matrix tform;//the current t-form, must be updated when something or parent changes

  Mat4x4 irest;//!<Inverse of the rest transform (armature space)
  Mat4x4  rest;//!<The rest transformation (armature space)

  Mat4x4 lmat; //!<The local matrix (loc/orient w.r.t. parent).
  Mat4x4 tform;//!<Cached complete animation tform (local to world)
  Mat4x4 lanim;//!<Complete local t-form: lmat*animation*length
  vector<Constraints> cons_cache;//position of end effector in this bones coord sys
  Joint * joint;

  Bone(){
    init((Bone*)0,Matrix::eye(4,4));
  }

  Bone(Bone * parent,
       const Matrix & lmat,
       double len=0.8,
       const char * name=0){
    init(parent,lmat,len,name);
  }

  void init(Bone * parent,
	    const Matrix & lmat,
	    double len=0.8,
	    const char * name=0){
    this->len=len;
    this->stiff=1;
    this->parent=0;

    this->selected=NotSelected;

    joint=new Joint();//EulerJointXYZ();

    if(parent)parent->addChild(this);
    
    cons_depth=0;

    ///this->lmat.fromMatrix(lmat);
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++)
	this->lmat(i,j)=lmat(i,j);
    }

    this->name=name;
    this->id=-1;//unknown id

    lanim.to_eye();
    //NOTE: set the rest and inverse rest transforms
    updateRestTransform(false);
    updateTransform();
  }
  void setJoint(Joint * joint){
    if(this->joint)delete this->joint;
    this->joint=joint;
    updateTransform(true);
  }
  bool connected(){
    return fabs(lmat(0,3))<1e-5 &&
      fabs(lmat(1,3))<1e-5 &&
      fabs(lmat(2,3))<1e-5;
  }
  Vec3d getPoseLoc() const {
    return joint->getPoseLoc();
  }
  Quaternion getPoseQuat() const{
    return joint->getPoseQuat();
  }
  void setPoseLoc(const Vec3d & p){
    joint->setPoseLoc(p);
  }
  void setPoseQuat(const Quaternion & q){
    joint->setPoseQuat(q);
  }
  void scale(double sc=1.2,bool recurse=true){
    int kids_connected=0;
    for(int i=0;i<(int)kids.size();i++){
      kids_connected+=kids[i]->connected();
    }
    if(kids_connected==(int)kids.size()){
      len*=sc;
      updateRestTransform(true);
      updateTransform();
    }
    else {
      printf("%s kids are not connected\n",name.c_str());
      double oldlen=len;
      len*=sc;
      updateRestTransform(true);
      updateTransform();
      for(int i=0;i<(int)kids.size();i++){
	Vec3d old(kids[i]->lmat(0,3),
		  kids[i]->lmat(1,3)+oldlen,
		  kids[i]->lmat(2,3));
	old*=sc;
	old.y-=len;

	kids[i]->lmat(0,3)=old.x;
	kids[i]->lmat(1,3)=old.y;
	kids[i]->lmat(2,3)=old.z;

	kids[i]->updateRestTransform(true);
	kids[i]->updateTransform();
      }
    }
    if(recurse)
      for(int i=0;i<(int)kids.size();i++)kids[i]->scale(sc);
  }
  void deselectAll(){
    selected=NotSelected;
    
    for(int i=0;i<(int)kids.size();i++)kids[i]->deselectAll();
  }
  vector<Bone *> assignID(){
    std::list<Bone*> q;
    q.push_back(this);

    int base=0;
    vector<Bone *> boneByID;
    while(!q.empty()){
      Bone * bone=q.front();
      q.pop_front();
      bone->id=base++;

      boneByID.push_back(bone);

      std::for_each(bone->kids.begin(),bone->kids.end(),
		    my_bind(q,&std::list<Bone*>::push_back));
      //	    std::mem_fun_bind1(q,&std::list<Bone*>::push_back));
    }
    return boneByID;
  }
  void addChild(Bone * kid){
    kids.push_back(kid);
    kid->parent=this;
  }
  void removeCons(){
    cons.n=0;
    for(int k=0;k<5;k++)
      cons.has_cons[k]=false;
  }
  void setCons(const Vec4d & c,int index){
    cons.set(index, c);
  }
  void drawSelect(){
    glPushMatrix();
    lmat.glMult();

    joint->glTransform();

    GLUquadric * q=gluNewQuadric();

    glPushName(id);

    glPushName(HeadSelected);//Head
    glPushMatrix();
    glColor4f(1,1,1,0.4);
    gluSphere(q,len/10,16,4);
    glPopMatrix();
    
    glLoadName(TipSelected);//Tail
    glPushMatrix();
    glTranslatef(0,len,0);
    glColor4f(1,1,1,0.4);
    gluSphere(q,len/14,16,4);
    glPopMatrix();
    
    /*
    glLoadName(BoneSelected);//Bone
    glBegin(GL_LINES);
    glVertex3f(0,len/10,0);
    glVertex3f(0,len-len/14,0);
    glEnd();
    */
    gluDeleteQuadric(q);
    glPopName();
    glPopName();
    
    glTranslatef(0,len,0);
    for(int i=0;i<(int)kids.size();i++){
      kids[i]->drawSelect();
    }
    glPopMatrix();
  }
  void draw(){
    glPushMatrix();
    //Matrix lmatt=lmat.transpose();
    //glMultMatrixd(lmatt.data);
    lmat.glMult();

    //glTranslatef(pos.x,pos.y,pos.z);
    //glRotatef(theta*180.0/M_PI,0,0,1);

    
    //Vec3f mx(0,cos(this->maxx)*len,sin(this->maxx)*len);

    /*
    int nlut=32;
    double * rlut=new double[nlut];
    for(int i=0;i<nlut;i++)rlut[i]=1;

    double miny=1;
    for(int i1=0;i1<=100;i1++){
      double p1=((double)i1/100.0)*(this->maxz-this->minz)+this->minz;
      
      for(int j1=0;j1<=100;j1++){
	double p2=((double)j1/100.0)*(this->maxx-this->minx)+this->minx;
	
	double x=-cos(p2)*sin(p1);
	double y= cos(p2)*cos(p1);
	double z= sin(p2);

	double theta=(atan2(x,z));
	int ti=int((theta+M_PI)/(2*M_PI)*(nlut));
	
	rlut[ti]=std::min(rlut[ti],y);
      }
    }
    glBegin(GL_LINE_LOOP);
    for(int i=0;i<nlut;i++){
      double phi=acos(rlut[i]);
      double theta=2.0*M_PI*(double)i/double(nlut);

      double x=sin(phi)*sin(theta)*len;
      double y=cos(phi)*len;
      double z=sin(phi)*cos(theta)*len;
      printf("y is %f %f\n",cos(phi),rlut[i]);
      glVertex3f(x,y,z);
    }
    glEnd();

    delete [] rlut;

    glRotatef(this->rz*180.0/M_PI,0,0,1);
    glRotatef(this->rx*180.0/M_PI,1,0,0);
    glRotatef(this->ry*180.0/M_PI,0,1,0);
    */
    if(Bone::drawBounds){
      //joint->drawBounds();
      /***
      bool single=(int(usex)+int(usey)+int(usez))==1;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      glDepthMask(0);

      glBegin(GL_LINE_STRIP);
      if(usex)
	for(double d=minx;d<=maxx;d+=std::max((maxx-minx)/20,0.01))
	  glVertex3f(0,cos(d)*len/2,sin(d)*len/2);
      if(usez && single)
	for(double d=minz;d<=maxz;d+=std::max((maxz-minz)/20,0.01))
	  glVertex3f(cos(d)*len/2,0,sin(d)*len/2);
      //no z as it is just roll
      glEnd();
      
      if(usex){
	glBegin(GL_TRIANGLE_FAN);
	glColor4f(0,0.9,0.9,0.5);
	glVertex3f(0,0,0);
	for(double d=minx;d<=rx;d+=std::max((rx-minx)/20,0.01)){
	  glVertex3f(0,cos(d)*len/2,sin(d)*len/2);
	}
	glEnd();
	
	glBegin(GL_TRIANGLE_FAN);
	glColor4f(0,0.4,0.4,0.5);
	glVertex3f(0,0,0);
	for(double d=rx;d<=maxx;d+=std::max((maxx-rx)/20,0.01)){
	  glVertex3f(0,cos(d)*len/2,sin(d)*len/2);
	}
	glEnd();
      }
      if(usez && single){
	glBegin(GL_TRIANGLE_FAN);
	glColor4f(0,0.9,0.9,0.5);
	glVertex3f(0,0,0);
	for(double d=minz;d<=rz;d+=std::max((rz-minz)/20,0.01)){
	  glVertex3f(cos(d)*len/2,0,sin(d)*len/2);
	}
	glEnd();
	
	glBegin(GL_TRIANGLE_FAN);
	glColor4f(0,0.4,0.4,0.5);
	glVertex3f(0,0,0);
	for(double d=rz;d<=maxz;d+=std::max((maxz-rz)/20,0.01)){
	  glVertex3f(cos(d)*len/2,0,sin(d)*len/2);
	}
	glEnd();
      }
      if(!single){
	glBegin(GL_LINE_LOOP);
	for(double d=this->miny;d<=this->maxy;d+=(this->maxy-this->miny)/100.0+0.0001){
	  double ct=cos(d);
	  double st=sin(d);

	  double xx=sin(this->maxx)*len*st;
	  double xy=cos(this->maxx)*len;
	  double xz=sin(this->maxx)*len*ct;

	  glVertex3f(xx,xy,xz);
	  
	}
	glEnd();
      }
      glDepthMask(1);
      ***/
    }
    joint->glTransform();

    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3f(0,0,0);
    glVertex3f(len/3,0,0);
    
    glColor3f(0,0,1);
    glVertex3f(0,0,0);
    glVertex3f(0,0,len/3);
    glEnd();

    glPushMatrix();

    glLineWidth(3);

    glColor3f(1,1,selected==BoneSelected?0.0:1.0);
    glBegin(GL_LINES);
    //glVertex3f(0,0,0);
    //glVertex3f(0,len,0);
    
    for(int k=0;k<2;k++){
      double xx=(k==0)?(len/8):0;
      double zz=(k==1)?(len/8):0;

      glColor4f(0.6,0.6,selected==HeadSelected?0.0:1.0,0.6);
      glVertex3f( xx,len/8,zz);
      glColor4f(0.8,0.8,selected==TipSelected?0.0:1.0,0.5);
      glVertex3f(0,len,0);
      
      glColor4f(0.6,0.6,selected==HeadSelected?0.0:1.0,0.6);
      glVertex3f(-xx,len/8,-zz);
      glColor4f(0.8,0.8,selected==TipSelected?0.0:1.0,0.5);
      glVertex3f(0,len,0);
      
      glColor4f(0.6,0.6,selected==HeadSelected?0.0:1.0,0.6);
      glVertex3f(0,0,0);
      glVertex3f(-xx,len/8,-zz);
      
      glColor4f(0.6,0.6,selected==HeadSelected?0.0:1.0,0.6);
      glVertex3f(0,0,0);
      glVertex3f( xx,len/8, zz);
    }
    glEnd();
    glLineWidth(1);
    
    glPopMatrix();

    GLUquadric * q=gluNewQuadric();

    /*glPointSize(6);
    glColor3f(0,0,1);
    glBegin(GL_POINTS);
    glVertex2f(0,0);
    glVertex2f(0,len);
    glEnd();
    glPointSize(1);
    */
    
    glPushMatrix();
    glColor4f(1,1,selected&HeadSelected?0.0:1.0,0.2);
    gluSphere(q,len/10,16,4);
    glPopMatrix();
    
    glTranslatef(0,len,0);
    glColor4f(1,1,selected&TipSelected?0.0:1.0,0.2);
    gluSphere(q,len/14,16,4);

    gluDeleteQuadric(q);

    glLineWidth(2);
    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex2f(0,0);
    glVertex2f(len/10.0,0);

    glColor3f(0,1,0);
    glVertex2f(0,0);
    glVertex2f(0,len/10.0);
    glEnd();
    glLineWidth(1);

    for(int i=0;i<(int)kids.size();i++)
      kids[i]->draw();

    glPopMatrix();
  }
  void drawConstraints(){
    glPushMatrix();

    GLubyte stip[1024];
    memset(stip,0xAA,1024);
    glPolygonStipple(stip);
    glEnable(GL_LINE_STIPPLE);
    glEnable(GL_POLYGON_STIPPLE);
    if(cons.has_cons[Bone::TipLocCons]){
      glColor4f(0.5,0.5,0.5,0.4);
      drawArrow(getPosition(),cons[Bone::TipLocCons]);
    }
    if(cons.has_cons[Bone::HeadLocCons]){
      glColor4f(0.5,0.5,0.5,0.4);
      drawArrow(getHead(),cons[Bone::HeadLocCons]);
    }
    glDisable(GL_LINE_STIPPLE);
    glDisable(GL_POLYGON_STIPPLE);
    
    for(int i=0;i<(int)kids.size();i++)kids[i]->drawConstraints();
    glPopMatrix();
  }
  void clearConstraintCache(){
    cons_cache.clear();
  }
  void updateCacheParent(const Constraints & cons,int depth){
    cons_cache.push_back(cons);
    
    if(parent && depth-1!=0){
      Constraints cpar;
      cpar=cons;
      for(int k=0;k<5;k++)
	if(cpar.has_cons[k])cpar[k]=lanim*cons.get(k);

      parent->updateCacheParent(cpar,depth-1);
    }
  }
  void updateCacheTip(){
    assert(cons.n);

    //FIXME: use constraint
    Constraints thecons;
    memcpy(thecons.has_cons,cons.has_cons,sizeof(cons.has_cons));
    thecons.n=cons.n;

    if(cons.has_cons[0])thecons[0]=Vec4d(0,0,0,cons_scales[0]);
    if(cons.has_cons[1])thecons[1]=Vec4d(0,-len * cons_scales[1],0,cons_scales[1]);
    if(cons.has_cons[2])thecons[2]=Vec4d(cons_scales[2],0,0,0);
    if(cons.has_cons[3])thecons[3]=Vec4d(0,cons_scales[3],0,0);
    if(cons.has_cons[4])thecons[4]=Vec4d(0,0,cons_scales[4],0);

    cons_cache.push_back(thecons);

    if(parent && cons_depth-1!=0){
      //multiply the constraint poses by the local transform
      for(int k=0;k<5;k++)
	if(thecons.has_cons[k])
	  thecons[k]=lanim*thecons[k];

      parent->updateCacheParent(thecons,cons_depth-1);
    }
  }
  void getTips(vector<Bone*> & tips){
    if(!kids.size()){
      tips.push_back(this);
      return;
    }
    for(int i=0;i<(int)kids.size();i++)
      kids[i]->getTips(tips);
  }
  void updateRestTransform(bool recurse=false){
    if(parent)rest=parent->rest;
    else rest.to_eye();
    
    rest*=lmat;
    rest.trans_right(0,len,0);
    irest=rest.eucInv();

    if(recurse)
      for(int k=0;k<(int)kids.size();k++)
	kids[k]->updateRestTransform(recurse);
  }
  void updateTransform(bool recurse=false){
    if(parent)tform=parent->getTransform();
    else tform.to_eye();

    lanim.to_eye();
    lanim.mult_right(lmat);//.copy();

    joint->apply(lanim);

    lanim.trans_right(0,len,0);

    tform.mult_right(lanim);

    if(recurse){
      for(int i=0;i<(int)kids.size();i++){
	kids[i]->updateTransform(recurse);
      }
    }
    
    updated++;
    //tform=tform*Matrix::trans(0,len,0);
    //if(parent)tform=parent->getTransform()*tform;//assume parent has been updated
  }
  void getDeriv(int pindex,int cindex,double * f,bool differential=false) const{
    Mat4x4 mat;
    if(parent){
      mat=parent->getTransform();
      mat.mult_right(lmat);
    }
    else mat=lmat;

    
    //int size=getNParams();//int(usex)+int(usey)+int(usez);
    
    Constraints dup=cons_cache[cindex];
    
    dup.trans(0,len,0);
    
    joint->applyDeriv(mat,pindex,differential);

    dup.mult_left(mat);

    int offs=0;
    for(int k=0;k<5;k++){
      //pt=mat*pt;
      if(dup.has_cons[k]){
	const Vec4d & pt=dup[k];
	f[offs+0]=pt.x;
	f[offs+1]=pt.y;
	f[offs+2]=pt.z;
	offs+=3;
      }
    }
  }
  void setDeriv(int index,bool differential=false){
    if(parent)tform=parent->getTransform();
    else tform.to_eye();
    tform.mult_right(lmat);//.copy();
    
    joint->applyDeriv(tform,index,differential);

    tform.trans_right(0,len,0);
  }
  Mat4x4 getRestTransform() const {
    return rest;
  }
  Mat4x4 getTransform() const{
    return tform;
  }
  Vec4d getPosition(){
    return Vec4d(tform(0,3),tform(1,3),tform(2,3),1);
  }
  Vec4d getHead(){
    return Vec4d(tform(0,3)-len*tform(0,1),
		 tform(1,3)-len*tform(1,1),
		 tform(2,3)-len*tform(2,1),1);
  }
  Vec4d getRotX(){
    return Vec4d(tform(0,0),tform(1,0),tform(2,0),0);
  }
  Vec4d getRotY(){
    return Vec4d(tform(0,1),tform(1,1),tform(2,1),0);
  }
  Vec4d getRotZ(){
    return Vec4d(tform(0,2),tform(1,2),tform(2,2),0);
  }
  int getNCons() const{
    if(cons.n)return 3*cons.n;
    return 0;
  }
  int getCons(double * f,bool deriv=false){
    assert(cons.n);

    int offs=0;
    Vec4d (Bone::*funcs [5])()={&Bone::getPosition,
				&Bone::getHead,
				&Bone::getRotX,
				&Bone::getRotY,
				&Bone::getRotZ};
    for(int k=0;k<5;k++){
      if(cons.has_cons[k]){
	Vec4d p=(this->*(funcs[k]))();
	f[offs+0]=(p.x-(deriv?0:cons[k].x)) * cons_scales[k];
	f[offs+1]=(p.y-(deriv?0:cons[k].y)) * cons_scales[k];
	f[offs+2]=(p.z-(deriv?0:cons[k].z)) * cons_scales[k];
	offs+=3;
      }
    }
    return offs;
  }
  int getNParamsAll(){
    int n=getNParams();
    for(int i=0;i<(int)kids.size();i++)n+=kids[i]->getNParamsAll();
    return n;
  }
  int getParamsFromAll(int * from,int * index){
    int offs=getNParams();
    for(int i=0;i<offs;i++){
      from[i]=id;
      index[i]=i;
    }
    for(int k=0;k<(int)kids.size();k++)offs+=kids[k]->getParamsFromAll(from+offs,index+offs);
    return offs;
  }
  int getParamsAll(double * p,bool differential=false){
    int offs=getParams(p,differential);
    for(int i=0;i<(int)kids.size();i++)
      offs+=kids[i]->getParamsAll(p+offs,differential);
    return offs;
  }
  int setParamsAll(const double * p,bool differential=false){
    int offs=setParams(p,differential,false);//NOTE: setParams has update parameter
    for(int i=0;i<(int)kids.size();i++)offs+=kids[i]->setParamsAll(p+offs,differential);
    return offs;
  }

  int getNParams() const{
    return joint->getNParams();
    ///return int(usex)+int(usey)+int(usez)+((useloc && !lockloc)?3:0);
  }

  int getParams(double * p,bool differential=false) const {
    return joint->getParams(p,differential);
  }

  int setParams(const double * p,bool differential,bool forceUpdate=false){
    int n=joint->setParams(p,differential);
    
    if(forceUpdate)updateTransform();//update before children (if something above has changed)
    return n;
  }

  vector<Bone *> getAllBones(){
    vector<Bone *> bones;
    bones.push_back(this);
    for(int i=0;i<(int)kids.size();i++){
      vector<Bone*> kidbones=kids[i]->getAllBones();
      bones.insert(bones.end(),kidbones.begin(),kidbones.end());
    }
    return bones;
  }
  /** \brief.  Get and return a list of inverse-kinematic chains.

      An inverse kinematic chain is a list of elements starting from
      the lowest link (with some constraint) and all the parents that affect
      this link.  This is done for each link with a constraint.
  */
  void getIKChains(vector<chain_t>& chains){
    static Timer funcTimer;
    static int called=0;
    if(!parent){
      funcTimer.reset();
      funcTimer.start();
    }
    if(cons.n){
      vector<Bone *> chain;
      Bone * bone=this;
      while(bone && ((int)chain.size()<cons_depth || !cons_depth)){
	chain.push_back(bone);
	bone=bone->parent;
      }
      chains.push_back(chain_t(constraint_t(this),chain));
    }
    for(int i=0;i<(int)kids.size();i++){
      kids[i]->getIKChains(chains);
    }
    //here we have all the IKChains
    if(!parent){
      called++;
     
      printf("there are %ld chains\n", chains.size());
      for(int i=0;i<(int)chains.size();i++){
	printf("%d: len=%ld\n", i, chains[i].second.size());
	for(int j=0;j<(int)chains[i].second.size();j++){
	  printf("\t%d %s\n", j, chains[i].second[j]->name.c_str());
	}
      }
     
      vector<Bone *> ikroots;//!< The root bones of the IK chains

      //Each IK root can be the root of several chains.  Collect these
      //and put all children in a list and keep track of all the constraints
      //and also keep track of what link affects which constraint.
      vector<vector<bone_affects_t> >   ikgroups;
      //For each group the ikaffects keeps track of constraints
      //of which constraints are part of a ikgroup
      vector<constraint_set_t> ikaffects;

      for(int i=0;i<(int)chains.size();i++){
	//The root is at the end of the chain
	Bone * theroot=chains[i].second[chains[i].second.size()-1];
	vector<Bone*>::iterator it=find(ikroots.begin(),ikroots.end(),theroot);

	//Try to see if we have something rooted at this bone
	if(it==ikroots.end()){
	  //In this case, we don't have a group rooted at this bone yet,
	  //so add one.
	  ikroots.push_back(theroot);

	  //All elements in the current chain affect ALL constraints in
	  //the tip.  Links are added parents first.
	  vector<bone_affects_t> chain;
	  for(int j=chains[i].second.size()-1;j>=0;j--){
	    chain.push_back(bone_affects_t(chains[i].second[j],0x1));
	  }
	  ikgroups.push_back(chain);
	  
	  //add the constraints to the list of constraints affected by this group
	  constraint_set_t cons_set;
	  cons_set.push_back(chains[i].first);
	  ikaffects.push_back(cons_set);
	}
	else{
	  //In this case, we already have a chain rooted at this group,
	  //so add the new constraints and any new links to the group

	  int ind=it-ikroots.begin();
	  //For every link in the chain, we check the group to see
	  //if the link already exists.  If it does, then we add
	  //a updated which constraints this link affects.  Otherwise,
	  //we have to add the element, and suggest that it only affects
	  //this latest constriant.
	  for(int j=chains[i].second.size()-1;j>=0;j--){
	    int found=-1;
	    for(int k=0;k<(int)ikgroups[ind].size();k++){
	      if(ikgroups[ind][k].first==chains[i].second[j]){
		found=k;
	      }
	    }
	    if(found>=0){
	      ikgroups[ind][found].second|=(0x1<<ikaffects[ind].size());
	    }
	    else{
	      ikgroups[ind].push_back(bone_affects_t(chains[i].second[j],(0x1<<ikaffects[ind].size())));
	    }
	  }

	  //assuming chain[0] doesn't already appear in the constraint set 
	  //for this group which it can't because bones are a tree.
	  ikaffects[ind].push_back(chains[i].first);
	}
      }
      
      printf("there are %ld roots\n",ikroots.size());
      for(int i=0;i<(int)ikroots.size();i++){
	printf("root %d at %s size %ld\n",i,ikroots[i]->name.c_str(),ikgroups[i].size());
	for(int j=0;j<(int)ikgroups[i].size();j++){
	  printf("\t%d %s affects: %X",j,ikgroups[i][j].first->name.c_str(),ikgroups[i][j].second);
	  printf("\n");
	}
      }
      bool differ=true;//false;//
      
      //should order the groups by depth from root first.
      //Now we can just solve the inverse kinematic problem for
      //each IK group.
      for(int gi=0;gi<(int)ikgroups.size();gi++){
	int m=0;//The number of constraints
	int n=0;//The number of degrees of freedom
	
	vector<int> cons_offsets;//Location of constraint in f-vector (0<loc<m)
	
	for(int j=0;j<(int)ikaffects[gi].size();j++){
	  cons_offsets.push_back(m);
	  m+=ikaffects[gi][j]->getNCons();
	}
	//For each element in the x-vector, keep track of which bone in the
	//group it came from and which parameter in that bone (index)
	vector<int> from; 
	vector<int> index;
	for(int j=0;j<(int)ikgroups[gi].size();j++){
	  int len=ikgroups[gi][j].first->getNParams();
	  if(len==0){
	    printf("%s has zero params\n",ikgroups[gi][j].first->name.c_str());
	  }
	  //keep track of which bone index this parameter came from
	  for(int k=0;k<len;k++){
	    from.push_back(j);
	    index.push_back(k);
	  }
	  n+=len;
	}
	printf("root: %s eq: %dx%d\n",ikroots[gi]->name.c_str(),m,n);

	Matrix fback(m,1);
	Matrix f(m,1);
	Matrix p(n,1);
	Matrix pback(n,1);
	
	//double * mins=new double[n];
	//double * maxs=new double[n];

	//Extract the parameters for this ik group, and the 
	//joint limits
	getParams(ikgroups[gi],p.data,differ);
	
	StopWatch itTimer;

	double t1 = 0;
	double solv = 0;

	for(int its=0;its<20;its++){
	  updated=0;

	  /***  
	  //project the rotation values on their joint limits
	  for(int i=0;i<n;i++){
	    p.data[i]=std::max(mins[i],p.data[i]);
	    p.data[i]=std::min(maxs[i],p.data[i]);
	  }
	  //force the use of the clamped parameters
	  setParams(ikgroups[gi],differ,p.data);
	  ***/
	  
	  //Backup parameters and get new ones (may have been clamped)
	  pback=p.copy();
	  getParams(ikgroups[gi],p.data,differ);
	  for(int i=0;i<n;i++)
	    if(pback[i]!=p[i])
	      printf("parameter %d for %s is clamped: %f,%f  %f\n",index[i],ikgroups[gi][from[i]].first->name.c_str(),p[i],pback[i],fabs(pback[i]-p[i]));

	  pback=p.copy();
	  
	  //get the values of all the constraints
	  getCons(ikaffects[gi],fback.data);
	  
	  
	  double res=sqrt(fback.dot(fback));
	  if(res<1e-6)break;
	  //printf("its %d res:%f\n",its,sqrt(fback.dot(fback)));

	  /*
	    Compute the Jacobian of the function.  We first
	    update the cached values of the constraints in the local
	    coordinate system of each bone.  Then just call getDeriv
	    on each element in the group to get the derivative of
	    the constraint w.r.t. to that bones parameters.
	  */
	  Matrix J(m,n);
	  J.setAll(0);
	  //double delta=1e-5;

	  Matrix J2(m,n);
	  J2.setAll(0);

	  StopWatch swatch;

	  for(int k=0;k<(int)ikgroups[gi].size();k++){
	    ikgroups[gi][k].first->clearConstraintCache();
	  }
	  //NOTE: used to clear here on first iteration,
	  //      but I think that would have caused bugs in some cases
	  for(int k=0;k<(int)ikaffects[gi].size();k++)
	    ikaffects[gi][k]->updateCacheTip();

	  //Now actually compute the jacobian
	  for(int i=0;i<n;i++){
	    //Only consider those parameters whose values are within
	    //the joint limits.
	    //if(pback[i]<mins[i]|| pback[i]>maxs[i])continue;
	    
	    uint32_t aff=ikgroups[gi][from[i]].second;
	    memset(f.data,0,sizeof(double)*m);
	    int which=0;//Used to index the constraint cache of the bone
	    for(int k=0;k<(int)ikaffects[gi].size();k++)
	      if(aff&(1<<k)){
		ikgroups[gi][from[i]].first->getDeriv(index[i],which,f.data+cons_offsets[k],differ);
		which++;
	      }

	    for(int j=0;j<m;j++)
	      J(j,i)=f[j];
	  }
	  t1+=double(swatch);
	  
	  swatch.reset();

	  /*
	  for(int i=0;i<n;i++){
	    //if(pback[i]<mins[i]||pback[i]>maxs[i])continue;
	    ikgroups[gi][from[i]].first->setDeriv(index[i]);

	    uint32_t aff=ikgroups[gi][from[i]].second;
	    //only update the transformation of those that are 
	    //below in the hierarchy and affect the constraint 
	    for(int k=from[i]+1;k<ikgroups[gi].size();k++)
	      if(ikgroups[gi][k].second & aff)
		ikgroups[gi][k].first->updateTransform();
	    		 
	    for(int k=0;k<ikaffects[gi].size();k++)
	      if(aff&(1<<k))
		ikaffects[gi][k].first->getCons(f.data+cons_offsets[k],true);

	    for(int j=0;j<m;j++)
	      J2(j,i)=f[j];
	    
	    //reset the transformations of this from[i] and its kids
	    for(int k=from[i];k<ikgroups[gi].size();k++)
	      if(ikgroups[gi][k].second & aff)
		ikgroups[gi][k].first->updateTransform();
	  }
	  t2+=double(swatch);
	  */
	  
	  /*
	  Matrix pabs=p.copy();
	  getParams(ikgroups[gi],pabs.data,false);//get absolute parameters

	  for(int i=0;i<n;i++){
	    //if(pback[i]<mins[i]||pback[i]>maxs[i])continue;
	    p[i]=pback[i]+delta;

	    setParams(ikgroups[gi],p.data,differ,from[i]);//only updates children
	
	    getCons(ikaffects[gi],f.data);
	    for(int j=0;j<m;j++)
	      J2(j,i)=(f[j]-fback[j])/delta;
	    
	    p[i]=pback[i];
	    setParams(ikgroups[gi],pabs.data,false,from[i],from[i]+1);
	    //setParams(ikgroups[gi],p.data,differ,from[i],from[i]+1);//only updates children
	  }
	  */
	  /*
	  J2.printMatlab("J2");
	  J.printMatlab("J");
	  J=J2;
	  printf("%f %f\n",t1,t2);
	  */
	  //J.printMatlab("J");
	  swatch.reset();

	  Matrix U,S,V;
	  J.Lsvd(U,S,V);
	  Matrix Sinv(V.m,U.m);
	  Sinv.setAll(0);
	  double alpha=0.5;
	  for(int i=0;i<std::min(Sinv.m,Sinv.n);i++){
	    Sinv(i,i)=(S[i]/(S[i]*S[i]+alpha*alpha));
	  }
	  fprintf(stderr,"Cond(%d): %f\n",its,S[0]/S[std::min(Sinv.m,Sinv.n)-1]);
	  Matrix VSinv(V.m,Sinv.n);
	  Matrix Jinv(U.m,V.m);
	  dgemm_('T','T',V.m,Sinv.n,V.n,1.0,V.data,V.n,Sinv.data,Sinv.n,0,VSinv.data,VSinv.m);
	  dgemm_('N','N',VSinv.m,U.m,U.n,1.0,VSinv.data,VSinv.m,U.data,U.n,0,Jinv.data,Jinv.n);
	  Jinv=Jinv.transpose();
	  //Matrix Jinv2=V*Sinv*U.transpose();
	  //Matrix diff=(Jinv-Jinv2);
	  //printf("norm:%f\n",diff.dot(diff));
	  
	  double len=sqrt(fback.dot(fback));
	  if(len>1)fback*=1/len;
  
	  Matrix upd=Jinv*fback;
	  p=pback-upd;
	  setParams(ikgroups[gi],p.data,differ);

	  //printf("updated: %d\n",updated);
	  solv+=double(swatch);
	}
	double total=double(itTimer);

	printf("total:%f  t1:%f rem:%f perc:%f\n",total,t1,total-t1,100.0*t1/total);
	printf("         solv:%f rem:%f perc:%f\n",solv,total-solv,100.0*solv/total);
	//printf("\nncons %d,nparams=%d\n",m,n);

	/***
	for(int i=0;i<n;i++){
	  p.data[i]=std::max(mins[i],p.data[i]);
	  p.data[i]=std::min(maxs[i],p.data[i]);
	}
	setParams(ikgroups[gi],p.data,differ);
	delete [] mins;
	delete [] maxs;
	***/
	ikroots[gi]->updateTransform(true);
      }
      double total=double(funcTimer);
      printf("total_func %f  %d  %s\n",total,called,(total>0.1)?"big":"");
    }
  }
  /** \brief Get the constraint_set and put the results in f.
      
  \param cons the constraints to retrieve.
  \param deriv boolean value indicating whether the request is for derivatives
   */
  void getCons(const constraint_set_t & cons,double * f,bool deriv=false){
    int n=0;
    for(int i=0;i<(int)cons.size();i++){
      n+=cons[i]->getCons(f+n,deriv);
    }
  }
  /** \brief Set the parameters in the group.
      
      \param group the group to update parameters
      \param p the parameters to use
      \param differential  flag differential mode
      \param start where to start the update
      \param end where to end the update 
      The forcing of updates is only done for the elements after start.
      This operation could be expensive if all the transformations always
      need to be updated.
   */
  int setParams(const vector<bone_affects_t> & group,
		double * p,
		bool differential,
		int start=0,int end=-1){
    int n=0;
    for(int i=0;i<(int)group.size();i++){
      if(end>=0 && i>=end)break;
      if(i>=start)n+=group[i].first->setParams(p+n,differential,i>=start);
      else n+=group[i].first->getNParams();
    }
    return n;
  }
  /** \brief Get the parameters for the group.

      \param group the group of parameters to get parameters from
      \param p the returned values of the parameters
      \params differential flag differential mode
  */
  int getParams(const vector<bone_affects_t> & group,
		double * p,bool differential){
    int n=0;
    for(int i=0;i<(int)group.size();i++)
      n+=group[i].first->getParams(p+n,differential);
    return n;
  }
  /** \brief Optimize the lengths of bones (these are constrained).
      
   */
  void optimizeLengths(vector<Bone*> &  boneByID,bool overall=false){
    std::list<Bone *> q;
    q.push_back(this);

    std::vector<Vec4d>  p;
    std::vector<Vec4d> dp;
    
    vector<int>  ids;
    vector<int>  inds;

    while(!q.empty()){
      Bone * bone=q.front();
      q.pop_front();

      if(bone->cons.has_cons[0]){
	p.push_back(bone->getPosition());
	dp.push_back(bone->cons[0]);

	ids.push_back(bone->id);
	inds.push_back(0);
      }
      if(bone->cons.has_cons[1]){
	p.push_back(bone->getHead());
	dp.push_back(bone->cons[1]);

	ids.push_back(bone->id);
	inds.push_back(1);
      }
      for(int i=0;i<(int)bone->kids.size();i++)
	q.push_back(bone->kids[i]);
    }
    if(!p.size())return;

    if(overall){
      Matrix A(3*p.size(),1);
      Matrix b(3*p.size(),1);
      for(int i=0;i<(int)p.size();i++){
	for(int j=0;j<3;j++){
	  A(3*i+j,0)=p[i].data[j];
	  b[3*i+j]=dp[i].data[j];
	}
      }
      Matrix sc=Matrix::LlinLeastSq(A,b);
      
      scale(sc[0],true);
    }
    //scale groups
    vector<vector<string> > scale_groups_s;
    vector<string> group;
    group.push_back("ltibia");
    group.push_back("rtibia");
    scale_groups_s.push_back(group);

    group.clear();
    group.push_back("lhumerus");
    group.push_back("rhumerus");
    scale_groups_s.push_back(group);
    
    group.clear();
    group.push_back("head");
    group.push_back("upperneck");
    scale_groups_s.push_back(group);

    group.clear();
    group.push_back("lowerback");
    group.push_back("upperback");
    group.push_back("lowerneck");
    group.push_back("thorax");
    scale_groups_s.push_back(group);

    group.clear();
    group.push_back("lhipjoint");
    group.push_back("rhipjoint");
    scale_groups_s.push_back(group);

    group.clear();
    group.push_back("lclavicle");
    group.push_back("rclavicle");
    scale_groups_s.push_back(group);

    group.clear();
    group.push_back("lradius");
    group.push_back("rradius");
    scale_groups_s.push_back(group);

    
    group.clear();
    group.push_back("lfemur");
    group.push_back("rfemur");
    scale_groups_s.push_back(group);

    vector<vector<int> > scale_groups;
    for(int i=0;i<(int)scale_groups_s.size();i++){
      vector<int> group;
      for(int j=0;j<(int)scale_groups_s[i].size();j++){
	int fnd=-1;
	for(int k=0;k<(int)boneByID.size();k++){
	  if(boneByID[k]->name==scale_groups_s[i][j]){
	    fnd=k;
	    break;
	  }
	}
	if(fnd>=0)group.push_back(fnd);
      }
      scale_groups.push_back(group);
    }

    bool * boneAffects=new bool[boneByID.size()*ids.size()];
    memset(boneAffects,0,boneByID.size()*ids.size());
    for(int i=0;i<(int)ids.size();i++){
      Bone * b=boneByID[ids[i]];
      while(b){
	boneAffects[b->id*ids.size()+i]=1;
	printf("%s affects %d\n",b->name.c_str(),i);
	b=b->parent;
      }
    }
    vector<vector<int> > groupAffects;

    vector<int> use_groups;
    for(int j=0;j<(int)scale_groups.size();j++){
      bool * merged=new bool[ids.size()];
      memset(merged,0,ids.size());
      for(int i=0;i<(int)scale_groups[j].size();i++){
	for(int k=0;k<(int)ids.size();k++)
	  merged[k]|=boneAffects[scale_groups[j][i]*ids.size()+k];
      }
      vector<int> affects;
      printf("%d affects: %ld/%ld: ",j,affects.size(),ids.size());
      for(int i=0;i<(int)ids.size();i++){
	if(merged[i]){
	  affects.push_back(i);
	  printf("%d ",i);
	}
      }
      groupAffects.push_back(affects);
      delete [] merged;
      printf("\n");
      
      if(affects.size())use_groups.push_back(j);
    }

    for(int its=0;its<10 && use_groups.size();its++){

      
      Matrix J(3*p.size(),use_groups.size());
      Matrix r(3*p.size(),1);
      for(int i=0;i<(int)p.size();i++){
	if(inds[i]==0)p[i]=boneByID[ids[i]]->getPosition();
	else p[i]=boneByID[ids[i]]->getHead();

	r[3*i  ]=(p[i].x-dp[i].x);
	r[3*i+1]=(p[i].y-dp[i].y);
	r[3*i+2]=(p[i].z-dp[i].z);
      }

      for(int j=0;j<(int)use_groups.size();j++){
	printf("use group #%d: %d\n",j,use_groups[j]);
	vector<int> & group=scale_groups[use_groups[j]];
	double eps=1.0/4096.0/4;
	double sc=1+eps;
	for(int k=0;k<(int)group.size();k++){
	  printf("%s \n",boneByID[group[k]]->name.c_str());
	  boneByID[group[k]]->scale(sc,0);
	}
	for(int k=0;k<(int)group.size();k++)
	  boneByID[group[k]]->updateTransform(true);
	updateTransform(true);

	Matrix f(J.m,1);
	f.setAll(0);
	for(int k=0;k<(int)groupAffects[use_groups[j]].size();k++){
	  int gi=groupAffects[use_groups[j]][k];
	  printf("affects %d\n",gi);
	  Vec4d pc;
	  if(inds[gi]==0)pc=boneByID[ids[gi]]->getPosition();
	  else pc=boneByID[ids[gi]]->getHead();

	  f[3*gi  ]=(pc.x-p[gi].x)/eps;
	  f[3*gi+1]=(pc.y-p[gi].y)/eps;
	  f[3*gi+2]=(pc.z-p[gi].z)/eps;
	}
	for(int k=0;k<f.m;k++)
	  J(k,j)=f[k];

	for(int k=0;k<(int)group.size();k++){
	  boneByID[group[k]]->scale(1.0/sc,0);
	}
	updateTransform(true);
      }
      Matrix upd=Matrix::LlinLeastSq(J,r);
      upd.printMatlab("upd");

      for(int j=0;j<(int)use_groups.size();j++){
	vector<int> & group=scale_groups[use_groups[j]];
	for(int k=0;k<(int)group.size();k++){
	  boneByID[group[k]]->scale(1.0-upd[j],0);
	}
	updateTransform(true);
      }
    }
    updateRestTransform(true);//update the rest transform if they changed

    vector<Bone::chain_t> chains;
    getIKChains(chains);
  }

  void setIKDepth(int depth,bool recurse=false){
    cons_depth=depth;
    if(recurse){
      depth=(depth==0)?0:(depth+1);
      for(int i=0;i<(int)kids.size();i++)
	kids[i]->setIKDepth(depth, true);
      //for_each(b->kids.begin(), b->kids.end(), boost::bind(Armature::setIKDepth,_1,(depth==0)?0:(depth+1), true));
    }
  }  
};

/*
bool operator==(const Bone::bone_affects_t & b1,
		const Bone::bone_affects_t & b2);
*/
#endif
