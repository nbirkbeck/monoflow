//#pragma interface

#ifndef ARMATURE_H
#define ARMATURE_H

#ifdef USING_GLEW
#include <GL/glew.h>
#endif

#include "bone.h"
#include "mesh.h"

#include <nmath/quaternion.h>
#include <nmath/vec3.h>
#include <vector>
#include <list>
#include <map>

using nacb::Quaternion;
using nacb::Vec3d;

struct JacobianCache{
  int nvert;
  int nparm;
  int nbone;
  int * index, * from;//!<length nparm
  bool * bone_affected;//!<boolean matrix (nbone x nparm) (row major) 
  vector<int> * affected_vert;

  JacobianCache(int nvert=0,int nparm=0,int nbone=0){
    bone_affected=0;
    index=from=0;
    affected_vert=0;
    alloc(nvert,nparm,nbone);
  }
  void alloc(int nvert,int nparm,int nbone){
    release();
    this->nvert=nvert;
    this->nparm=nparm;
    this->nbone=nbone;

    bone_affected=new bool[nbone*nparm];
    memset(bone_affected,0,sizeof(bool)*nbone*nparm);
    
    index=new int[nparm];
    from=new int[nparm];
    
    affected_vert=new vector<int>[nbone];//<This should be nbone, shouldn't it? Was nparm
  }
  void release(){
    nvert=nparm=nbone=0;

    delete [] bone_affected;
    delete [] index;
    delete [] from;
    delete [] affected_vert;

    bone_affected=0;
    affected_vert=0;
    index=from=0;
  }
  ~JacobianCache(){
    release();
  }
};

class Armature;
class PoseKeys{
 private:
  vector<vector<nacb::Quaternion> > quats;
  vector<vector<nacb::Vec3d> > pos;
  int current;
  int nframes;
  std::map<int, int> posIndex; //!< Map a global id to the restricted pos vecetor.

 public:
  vector<string> names;

  PoseKeys(){
    current=0;
    nframes=0;
  }

  PoseKeys(const vector<string> & names){
    this->names = names;
    current = 0; //FIXME:-1?
    nframes = 0;
  }

  PoseKeys resample(int skip) {
    PoseKeys keys = *this;
    
    for (int i=0; i<(int)quats.size(); i++) {
      keys.quats[i].clear();

      for (int k=0; k<(int)quats[i].size(); k+=skip) {
	keys.quats[i].push_back(quats[i][k]);
      }
    }
    
    for (int i=0; i<(int)pos.size(); i++) {
      keys.pos[i].clear();

      for (int k=0; k<(int)pos[i].size(); k+=skip) {
	keys.pos[i].push_back(pos[i][k]);
      }
    }    
    
    keys.nframes = keys.quats[0].size();
    return keys;
  }

  PoseKeys stripFront(int n) {
    PoseKeys keys = *this;
    
    printf("Removing: %d from quats.\n", n);
    for (int i=0; i<(int)quats.size(); i++) {
      keys.quats[i] = std::vector<nacb::Quaternion>(quats[i].begin() + n, quats[i].end());
    }
    
    printf("Removing: %d from pos.\n", n);
    for (int i=0; i<(int)pos.size(); i++) {
      if (pos[i].size())
	keys.pos[i] = std::vector<nacb::Vec3d>(pos[i].begin() + n, pos[i].end());
    }      
    
    keys.nframes -= n;
    return keys;
  }

  void append(Armature * arm);

  void append(const vector<nacb::Quaternion> & q,
	      const vector<nacb::Vec3d> & p,
	      const std::vector<int>& withPos = std::vector<int>()) {
    assert(!quats.size() || quats.size() == q.size());
    assert(!pos.size() || pos.size() == p.size());

    printf("%ld %ld\n", posIndex.size(), withPos.size());
    assert(!posIndex.size() || posIndex.size() == withPos.size());

    if (!posIndex.size()) {
      for (int i=0; i<(int)withPos.size(); i++) {
	posIndex[withPos[i]] = i;
      }
    }
    else {
      for (int i=0; i<(int)withPos.size(); i++) {
	if (!posIndex.count(withPos[i])) {
	  throw std::runtime_error("Mapping has changed (key doesnt exist)");
	}
	else if (posIndex[withPos[i]] != i) {
	  throw std::runtime_error("Mapping has changed!");
	}
      }
    }

    if(!quats.size())
      quats = vector<vector<nacb::Quaternion> >(q.size());
    if(!pos.size())
      pos = vector<vector<nacb::Vec3d> >(p.size());

    for(int i=0; i<(int)quats.size(); i++){
      quats[i].push_back(q[i]);
    }
    for(int i=0; i<(int)pos.size(); i++){
      pos[i].push_back(p[i]);
    }
    nframes = quats[0].size();
  }

  int save(const std::string & str) const {
    return save(str.c_str());
  }

  int save(const char * filename) const {
    FILE * file = fopen(filename, "w");
    if(!file)return 0;

    fprintf(file, "ANIM %d\n", 1);

    int nbones = names.size();
    fprintf(file, "%d\n", nbones);
    
    for(int i=0; i<(int)names.size(); i++)
      fprintf(file, "%s\n", names[i].c_str());
    
    fprintf(file, "%d\n", nframes);

    // Begin (Version > 0)
    {
      std::map<int, int>::const_iterator it(posIndex.begin());
      std::vector<std::pair<int, int> > posIndexVector;
      
      // Map is sorted by key, we need to sort by value (e.g., the index)
      // Create a list of element pairs (key is the value of the map).
      for (; it != posIndex.end(); ++it) {
	posIndexVector.push_back(make_pair(it->second, it->first));
      }
      std::sort(posIndexVector.begin(), posIndexVector.end());
      
      
      fprintf(file, "%ld ", posIndex.size());
      for (int i=0; i<(int)posIndexVector.size(); i++) {
	fprintf(file, "%d ", posIndexVector[i].second);
      }
      fprintf(file, "\n");
    }
    // End (version > 0)

    assert(pos.size()>0);

    for (int i=0; i<(int)quats.size(); i++) {
      // Output the pos first (if we have it).
      int pi = posIndex.count(i) ? posIndex.find(i)->second : -1;
      if (pi >= 0) {
	for(int i=0; i<nframes; i++){
	  fprintf(file, "%g %g %g ", pos[pi][i].x, pos[pi][i].y, pos[pi][i].z);
	}
	fprintf(file, "\n");
      }

      // Output the quaternion.
      for(int j=0; j<nframes; j++){
	fprintf(file, "%g %g %g %g ", quats[i][j].v.x, 
		quats[i][j].v.y, quats[i][j].v.z, quats[i][j].a);
      }
      fprintf(file, "\n");
    }
    fclose(file);
    return nframes;
  }

  // Convenience function.
  int load(const std::string & str){
    return load(str.c_str());
  }

  int load(const char * filename){
    FILE * file=fopen(filename,"r");
    if(!file)return 0;
    int nbones;
    char line[1024] = "";

    posIndex.clear();
    pos.clear();
    quats.clear();

    int version = 0; 
    //Strip the comments from the beginning of the file.
    while(fgets(line, 1024, file) && line[0] == '#'){
      printf("Stripping comment: %s\n", line);
    }
    
    if (strstr(line, "ANIM") == line) {
      if (1 != sscanf(line + 4, "%d", &version)) {
	printf("Error loading version.\n");
      }
      printf("Got ANIM version: %d\n", version);
      fgets(line, 1024, file);
    }

    sscanf(line,"%d",&nbones);

    for(int i=0;i<nbones;i++){
      fgets(line,1024,file);
      if(strlen(line) && line[strlen(line)-1]=='\n')
	line[strlen(line)-1]=0;
      names.push_back(line);

      quats.push_back(vector<nacb::Quaternion>());

      printf("PoseKeys: %s\n",line);
    }
    fscanf(file,"%d",&nframes);

    if (version == 0) {
      // Just one position (for root).
      pos.push_back(vector<nacb::Vec3d>());
      posIndex[0] = 0;
    }
    else {
      int numWithPos = 0;
      fscanf(file, "%d", &numWithPos);

      if (numWithPos < 0 || numWithPos > nbones) {
	fprintf(stderr, "There are no position keys (or there are too many)\n");
	return 0;
      }
      
      // Add position keys for each element with a pose.
      for (int i=0; i<numWithPos; i++) {
	int index;
	fscanf(file, "%d", &index);
	pos.push_back(vector<nacb::Vec3d>());
	
	posIndex[index] = i;
      }
    }
    //printf("PoseKeys: Got LocKeys %d, nframes: %d\n",pos[0].size(), nframes);

 
    for(int i=0; i<nbones; i++){
      int pi = posIndex.count(i) ? posIndex[i] : -1;
      if (pi >= 0) {
	for(int i=0; i<nframes; i++){
	  Vec3d p;
	  if(3==fscanf(file,"%lf %lf %lf",&(p.x),&(p.y),&(p.z))){
	    pos[pi].push_back(p);
	    cout << "PoseKey: " << p << endl;
	  }
	}
	printf("PoseKeys: Got LocKeys %ld, nframes: %d\n",pos[pi].size(), nframes);
      }
      nacb::Vec3d v;
      double a;
      for(int j=0;j<nframes;j++){
	if(4==fscanf(file,"%lf %lf %lf %lf",&(v.x),&(v.y),&(v.z),&a))
	  quats[i].push_back(nacb::Quaternion(v,a));
      }
      //printf("PoseKeys: got %d keys for bone %d\n",quats[i].size(),i);
    }

    //printf("PoseKeys: read in keys for %d bones and nframes=%d\n",nbones,nframes);

    current=0;

    fclose(file);
    return nframes;
  }

  bool next(){
    current++;
    return (current<nframes);
  }

  Quaternion getQuat(int i){
    Quaternion q;
    if(i<0 || i>=(int)quats.size())return q;
    if(current<0 || current>=(int)quats[i].size())return q;
    return quats[i][current];
  }

  Vec3d getLoc(int i){
    Vec3d p(0,0,0);

    int posi = posIndex.count(i) ? posIndex[i] : -1;

    if(posi<0 || posi>=(int)pos.size())return p;
    if(current<0 || current>=(int)pos[posi].size())return p;

    //cout << "PoseKey: return " << pos[i][current] << endl;
    return pos[posi][current];
  }

  Vec3d setLoc(int i, const Vec3d & value){
    int posi = posIndex.count(i) ? posIndex[i]: - 1;
    if(posi >= (int)pos.size() || current >= (int)pos[posi].size()){
      throw "Pos index out of bounds";
    }
    pos[posi][current] = value;
    return pos[posi][current];
  }

  int getCurrent() const {
    return current;
  }

  int getNumFrames() const {
    return nframes;
  }
  
  bool setCurrent(int cur) {
    current = cur;
    return (current<nframes && current>=0);
  }
};


class Armature{
 protected:
  vector<Bone *> boneByID;
  std::map<std::string,Bone*> boneByName;

 public:
  vector<Bone *> tips; //!< List of elements that can be used for inverse kinematics
  Bone * root;
  int nparms;
  int *  from;
  int * index;

  Mesh * mesh;
  
  Armature(Bone * root=0){
    this->root = 0;
    from = index = 0;
    nparms = 0;

    if(root)init(root);
  }  
  ~Armature(){
    release();
  }
  void release(){
    root=0;//NOTE: doesn't delete root...ever
    boneByID.clear();
    boneByName.clear();
    tips.clear();

    delete [] from;
    delete [] index;
    from=0;
    index=0;
    nparms=0;
  }
  void init(Bone * root){
    assert(root);

    release();

    this->root=root;
    boneByID=root->assignID();
    for(int i=0;i<(int)boneByID.size();i++){
      boneByName[boneByID[i]->name]=boneByID[i];
    }
    if(root){
      nparms=root->getNParamsAll();
      from = new int[nparms];
      index = new int[nparms];
      
      root->getParamsFromAll(from,index);
    }
    else {
      nparms = 0;
      from = index = 0;
    }

    tips=boneByID;//same
  }
  
  void initBasicArm(int nsegs=5){
    Bone * parent=0;
    Bone * theRoot=0;
    //Bone * middle=0;
    int bonei=0;
    char name[1024];
    for(int i=0;i<nsegs;i++){
      sprintf(name,"Bone%d",bonei);
      parent=new Bone(parent,Matrix::eye(4,4),0.4,name);
      parent->setJoint(new EulerJointXYZ());
      if(!theRoot)theRoot=parent;
      //if(i==1)middle=parent;
      bonei++;
    }

    /*parent=middle;
    for(int i=0;i<nsegs/2;i++){
      sprintf(name,"Bone%d",bonei);
      parent=new Bone(parent,Matrix::eye(4,4),0.4,name);
      bonei++;
      }*/

    init(theRoot);
  }
  
  void useEuclideanRoot(){
    if (dynamic_cast<EulerJoint*>(root->joint)) {
      // We don't want to change the type of joint
      root->joint->useloc=true;
      root->joint->lockloc=false;
      setRootJoint(root->joint);
      return;
    }
    EulerJointXYZ * joint=new EulerJointXYZ(true,true,true);
    joint->useloc=true;
    joint->lockloc=false;
    setRootJoint(joint);
    printf("using euclidean root (nparms=%d, nbones=%ld)\n",
           nparms, boneByID.size());
  }

  void setRootJoint(Joint * joint) {
    if (joint != root->joint)
      root->setJoint(joint);

    // Need to re-allocate the from/index and nparms
    delete [] from;
    delete [] index;

    nparms = root->getNParamsAll();
    from = new int[nparms];
    index = new int[nparms];
  }

  void translate(int ww,int wh,double fov,Vec3d & cpos,Quaternion & q,double dx,double dy){
    double aspect=double(ww)/double(wh);
    Mat4x4 model(cpos,q);
    Vec3d xdir(model(0,0),model(0,1),model(0,2));
    Vec3d ydir(model(1,0),model(1,1),model(1,2));

    xdir*=2*aspect*dx*tan(fov*M_PI/180.0/2.0)*cpos.z/double(ww);
    ydir*=2*dy*tan(fov*M_PI/180.0/2.0)*cpos.z/double(wh);
    
    Vec3d t=ydir-xdir;
    Vec4d lt=(root->lmat.eucInv())*Vec4d(t.x,t.y,t.z,0.0);
    
    Vec3d was=root->getPoseLoc();
    printf("was %f %f %f\n",was.x,was.y,was.z);
    root->setPoseLoc(root->getPoseLoc()+Vec3d(lt.x,lt.y,lt.z));
    
    printf("loc changed %f %f %f\n",lt.x,lt.y,lt.z);
    
    root->updateTransform(true);
  }
  void rotate(int ww,int wh,double fov,Vec3d & cpos,Quaternion & qcam,double dx,double dy){
    Mat4x4 model(cpos,qcam);
    nacb::Vec3d axis(model(2,0),model(2,1),model(2,2));
    double angle=2*M_PI*dx/double(ww);      
    Mat4x4 mat=Mat4x4::angleAxis(axis,angle);
    
    Mat4x4 tform=root->lmat;
    root->joint->apply(tform);
    
    Vec3d tr(tform(0,3),tform(1,3),tform(2,3));
    tform(0,3)=tform(1,3)=tform(2,3)=0;
    
    mat=Mat4x4::trans(tr.x,tr.y,tr.z)*mat*tform;
    
    mat=root->lmat.eucInv()*mat;
    Quaternion q=Quaternion::fromMatrix(mat);
    root->setPoseQuat(q);
    root->setPoseLoc(Vec3d(mat(0,3),mat(1,3),mat(2,3)));
    
    root->updateTransform(true);
  }

  void optimizeLengths(){
    root->optimizeLengths(boneByID);
  }

  double guessScale(Mesh & mesh){
    double scale=1;

    if(mesh.hasGeometry()){
      double miny= 10e10;
      double maxy=-10e10;
      for(int i=0;i<(int)boneByID.size();i++){
	Vec4d p=boneByID[i]->getPosition();
	miny=std::min(miny,p.y);
	maxy=std::max(maxy,p.y);
	
	p=boneByID[i]->getHead();
	miny=std::min(miny,p.y);
	maxy=std::max(maxy,p.y);
      }
      double mesh_miny= 10e10;
      double mesh_maxy=-10e10;
      for(int i=0;i<(int)mesh.restVert.size();i++){
	mesh_miny=std::min((double)mesh.restVert[i].y,mesh_miny);
	mesh_maxy=std::max((double)mesh.restVert[i].y,mesh_maxy);
      }
      scale=0.9*(mesh_maxy-mesh_miny)/(maxy-miny);
      printf("scale %f\n",scale);
      ///scaledArmature=scale;
      root->scale(scale,true);
      
      //root->lmat(1,3)+=mesh_miny-miny*scale+(mesh_maxy-mesh_miny)*0.05;
      root->setPoseLoc(Vec3d(0,mesh_miny-miny*scale+(mesh_maxy-mesh_miny)*0.05,0));
      root->updateRestTransform(true);
      root->updateTransform(true);
    }
    else{
      double miny= 10e10;
      double maxy=-10e10;
      for(int i=0;i<(int)boneByID.size();i++){
	Vec4d p=boneByID[i]->getPosition();
	miny=std::min(miny,p.y);
	maxy=std::max(maxy,p.y);
	
	p=boneByID[i]->getHead();
	miny=std::min(miny,p.y);
	maxy=std::max(maxy,p.y);
      }
      scale=2.0/(maxy-miny+0.0001);
      printf("scale is %f\n",scale);
      
      root->scale(scale,true);
      root->setPoseLoc(Vec3d(0,-miny*scale,0));
      root->updateRestTransform(true);
      root->updateTransform(true);
    }
    return scale;
  }

  int read(const char * filename){
    release();
    FILE * file=fopen(filename,"r");
    if(!file)
      return 0;
    
    char line[2048];
    vector<Bone *> bones;

    Bone * theRoot=0;

    while(fgets(line,2048,file)){
      char name[256];
      Matrix mat=Matrix::eye(4,4);
      int id,parent;
      double len;int flip;
      
      /*int succ=*/
      sscanf(line,"%d %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d ",
             &id,name,&parent,mat.data+0,mat.data+1,mat.data+2,mat.data+3,mat.data+4,mat.data+5,mat.data+6,mat.data+7,mat.data+8,mat.data+9,mat.data+10,mat.data+11,&len,&flip);
      
      while(id>=(int)bones.size())bones.push_back(0);
      
      Bone * parentBone=(parent>=0)?bones[parent]:0;
      
      //Don't scale here anymore, it is a command line argument
      /*
      double scale=0.1;
      mat(0,3)*=scale;
      mat(1,3)*=scale;
      mat(2,3)*=scale;
      len*=scale;
      */
      char * mit=strstr(line, " mit ");
      char * vx=strstr(line, "vx(");
      char * vy=strstr(line, "vy(");
      char * vz=strstr(line, "vz(");

      char * rx=strstr(line,"rx(");
      char * ry=strstr(line,"ry(");
      char * rz=strstr(line,"rz(");
      char * aa=strstr(line, "aa()");
      char * pos=strstr(line, "pos()");
      float minx=0,maxx=0,miny=0,maxy=0,minz=0,maxz=0;

      std::vector<EulerJoint::Dof> dof;
      //The comma is actually what I wanted, but I accidentally exported
      //some files with a space, hence the second case!
      if(rx){
	if(2!=sscanf(rx,"rx(%f,%f)",&minx,&maxx))
	  sscanf(rx,"rx(%f %f)",&minx,&maxx);
      }
      if(ry){
	if(2!=sscanf(ry,"ry(%f,%f)",&miny,&maxy))
	  sscanf(ry,"ry(%f %f)",&miny,&maxy);
      }
      if(rz){
	if(2!=sscanf(rz,"rz(%f,%f)",&minz,&maxz))
	  sscanf(rz,"rz(%f %f)",&minz,&maxz);
      }
      
      // mit has a pose as well.
      double valx = 0, valy = 0, valz = 0;
      if(vx && 1 == sscanf(vx, "vx(%lf)", &valx))
	;
      if(vy && 1 == sscanf(vy, "vy(%lf)", &valy))
	;
      if(vz && 1 == sscanf(vz, "vz(%lf)", &valz))
	;

      const double d2r=M_PI/180.0;

      if(rx) 
	dof.push_back(EulerJoint::Dof(0, valx*d2r, minx*d2r, maxx*d2r));
      if(ry) 
	dof.push_back(EulerJoint::Dof(1, valy*d2r, miny*d2r, maxy*d2r));
      if(rz) 
	dof.push_back(EulerJoint::Dof(2, valz*d2r, minz*d2r, maxy*d2r));

      bones[id]=new Bone(parentBone,mat,len,name);
      
      Joint * joint = 0;
      if(aa)
	joint = new RodJoint();
      else if(mit)
	joint = new EulerJoint(dof);
      else
	joint = new EulerJointXYZ(minx<maxx,
				  miny<maxy,
				  minz<maxz,
				  minx*d2r,maxx*d2r,
				  miny*d2r,maxy*d2r,
				  minz*d2r,maxz*d2r);
      joint->useloc=false;
      joint->lockloc=true;
      if (pos != 0) {
	joint->useloc = true;
	joint->lockloc = false;
      }
      bones[id]->setJoint(joint);
      /***
      bones[id]->minx=minx*M_PI/180.0;
      bones[id]->maxx=maxx*M_PI/180.0;
      bones[id]->miny=miny*M_PI/180.0;
      bones[id]->maxy=maxy*M_PI/180.0;
      bones[id]->minz=minz*M_PI/180.0;
      bones[id]->maxz=maxz*M_PI/180.0;
      
      bones[id]->usex=minx<maxx;
      bones[id]->usey=miny<maxy;
      bones[id]->usez=minz<maxz;
      ***/
      printf("%f %f  %f %f  %f %f\n",minx,maxx,miny,maxy,minz,maxz);

      if(!parentBone)theRoot=bones[id];
    }
    
    fclose(file);

    init(theRoot);
    return true;
  }

  
  void setIKDepths(const std::vector<std::string> & names,
		   const std::vector<int> & depths,
		   bool recurse=false){
    assert(names.size()==depths.size());
    
    std::vector<std::string>::const_iterator b=names.begin();
    std::vector<int>::const_iterator d=depths.begin();
    
    while(b!=names.end()){
      Bone * bone=operator[](*b);
      if(bone)
	bone->setIKDepth(*d,recurse);
      else
	printf("cannot find bone with name %s\n",(*b).c_str());
      b++;
      d++;
    }
  }
  size_t size(){
    return boneByID.size();
  }

  const Bone * operator[](int id) const {
    if(id<(int)boneByID.size() && id>=0)
      return boneByID[id];
    return 0;
  }

  Bone * operator[](int id){
    if(id<(int)boneByID.size() && id>=0)
      return boneByID[id];
    return 0;
  }

  Bone * operator[](const std::string & name){
    return boneByName[name];
  }

  void set(PoseKeys & pk){
    for(int i=0;i<(int)pk.names.size();i++){
      Bone * bone=operator[](pk.names[i]);
      //printf("bone is %p (%s)\n",bone,pk.names[i].c_str());
      if(!bone)continue;
      bone->setPoseQuat(pk.getQuat(i));
      bone->setPoseLoc(pk.getLoc(i));

      //printf("setting bone %d\n",i);
    }
    root->updateTransform(true);
  }
  
  std::vector<Mat4x4> getTransforms() const{
    std::vector<Mat4x4>    tforms;
    
    root->updateTransform(true);
    
    for(int i=0;i<(int)boneByID.size();i++){
      tforms.push_back(boneByID[i]->tform*boneByID[i]->irest);
    }
    return tforms;
  }

  void animate(Mesh & mesh){
    root->updateTransform(true);
    
    Mat4x4 * tforms=new Mat4x4[boneByID.size()];
    for(int i=0;i<(int)boneByID.size();i++){
      tforms[i]=boneByID[i]->tform*boneByID[i]->irest;
    }
    mesh.animate(tforms);
    delete [] tforms;
  }

  void animate(Mesh & mesh,double * animp,int animp_len){
    double * back=new double[animp_len];
    root->getParamsAll(back);

    for(int i=0;i<animp_len;i++){
      animp[i]=(double(rand())/double(RAND_MAX))*0.25;
    }
    root->setParamsAll(animp);
    root->updateTransform(true);
    
    Mat4x4 * tforms=new Mat4x4[boneByID.size()];
    for(int i=0;i<(int)boneByID.size();i++){
      tforms[i]=boneByID[i]->tform*boneByID[i]->irest;
    }
    mesh.animate(tforms);
    
    root->setParamsAll(back);
    root->updateTransform();
    
    delete [] tforms;
    
    delete [] back;
  }


  //FIXME:
  void synchMesh(){

  }

  void relaxHips(bool test=false);

  void updateTransform(){
    root->updateTransform(true);
  }

  void setParams(const double * parms,bool differ=false){
    root->setParamsAll(parms,differ);
  }
  int getParams(double * parms,bool differ=false){
    return root->getParamsAll(parms,differ);
  }
  int getNParams() const{
    return nparms;
  }

  double * meshJacobian_fd(Mesh & mesh,
			   const double * animp,
			   int animp_len,bool differential=false);
  
  double * meshJacobian(Mesh & mesh,
			const double * animp,
			int animp_len,
			bool differential=false,
			JacobianCache * jout=0);

  void applyCollisionConstraint(Mesh & mesh);

  std::vector<std::string> getBoneNames() {
    std::vector<std::string> names;
    for(int i=0; i<(int)size(); i++){
      names.push_back(boneByID[i]->name);
    }
    return names;
  }

  bool writeArmature(const char * filename){
    FILE * file=fopen(filename,"w");
    if(!file)return false;
    
    for(int i=0;i<(int)boneByID.size();i++){
      Bone * b=boneByID[i];
      //FIXME: flip bool is wrong
      fprintf(file,"%d %s %d %g %g %g %g %g %g %g %g %g %g %g %g %g %d ",
	      i,b->name.c_str(),b->parent?b->parent->id:-1,
	      b->lmat(0,0),b->lmat(0,1),b->lmat(0,2),b->lmat(0,3),
	      b->lmat(1,0),b->lmat(1,1),b->lmat(1,2),b->lmat(1,3),
	      b->lmat(2,0),b->lmat(2,1),b->lmat(2,2),b->lmat(2,3),b->len,0);
      
      b->joint->writeBounds(file);
      
      /***	
		double rad2deg=180.0/M_PI;
		if(b->usex)fprintf(file,"rx(%g,%g) ",b->minx*rad2deg,b->maxx*rad2deg);
		if(b->usey)fprintf(file,"ry(%g,%g) ",b->miny*rad2deg,b->maxy*rad2deg);
		if(b->usez)fprintf(file,"rz(%g,%g) ",b->minz*rad2deg,b->maxz*rad2deg);
      ***/
      fprintf(file,"\n");
    }
    fclose(file);
    return true;
  }
  bool writePose(const char * filename){
    FILE * file=fopen(filename,"w");
    if(!file)return false;
    fprintf(file,"APOS\n");
    for(int i=0;i<(int)boneByID.size();i++){
      Quaternion q=boneByID[i]->getPoseQuat();
      fprintf(file,"%s QUAT %g %g %g %g\n",boneByID[i]->name.c_str(),q.v.x,q.v.y,q.v.z,q.a);
      //NOTE: remove this access
      if(boneByID[i]->joint->useloc){
	Vec3d l=boneByID[i]->getPoseLoc();
	fprintf(file,"%s POS %g %g %g\n",boneByID[i]->name.c_str(),l.x,l.y,l.z);
      }
    }
    fclose(file);
    return true;
  }  
  bool readPose(const char * filename){
    printf("reading pose from %s\n",filename);

    FILE * file=fopen(filename,"r");
    char line[2048],name[512];
    if(!fgets(line,1024,file))return false;
  
    //APOS
    if(strcmp(line,"APOS\n")){
      printf("file %s doesnt start with APOS\n",filename);
      fclose(file);
      return false;
    }
    while(fgets(line,2048,file)){
      double q[4];
      if(5==sscanf(line,"%s QUAT %lf %lf %lf %lf",name,q,q+1,q+2,q+3)){
	printf("%s quat %f %f %f %f\n",name,q[0],q[1],q[2],q[3]);
	Bone * bone=0;
	for(int i=0;i<(int)boneByID.size();i++){
	  if(boneByID[i]->name==name){
	    bone=boneByID[i];
	    break;
	  }
	}
	if(bone){
	  Quaternion qt=Quaternion(Vec3d(q[0],q[1],q[2]),q[3]);
	  qt.normalize();
	  bone->setPoseQuat(qt);
	}
      }
      else if(4==sscanf(line,"%s POS %lf %lf %lf",name,q+0,q+1,q+2)){
	//FIXME: wow this is sketchy as shit, and is actually wrong when root->lmat(0:3,0:3) is not eye(3,3)
	//assert(name==root->name);
	operator[](name)->setPoseLoc(Vec3d(q[0],q[1],q[2]));
      }
      else if(4==sscanf(line,"%s RXYZ %lf %lf %lf",name,q+0,q+1,q+2)){
	printf("%s rx:%f ry:%f rz:%f\n",name,q[0],q[1],q[2]);
      }
    }
    root->updateRestTransform(true);
    root->updateTransform(true);

    fclose(file);
    return true;
  }

  bool loadobj(const std::vector<std::string> & strings){
    release();

    std::vector<Bone *> boners;

    Bone * theRoot = 0;
    for(int i=0; i<(int)strings.size(); i++){
      if(strings[i].find("b ") == 0){
	std::stringstream ss(strings[i]);
       	std::string bn, name;
	int pid;
	Quaternion q;
	double len;
	Vec3d p;

	printf("got bone string: %s\n", strings[i].c_str());fflush(stdout);

	ss >> bn >> name >> pid >> len >> q.v.x >> q.v.y >> q.v.z >> q.a >> p.x >> p.y >> p.z;
	
	Bone * parent = 0;
	if(pid>=0 && pid<(int)boners.size())
	  parent = boners[pid];
	
	Mat4x4 lmat(p, q);
	Matrix mat(4,4);
	memcpy(mat.data, lmat.data, sizeof(double)*16);
	Bone * bone = new Bone(parent, mat, len, name.c_str());
	bone->setJoint(new EulerJointXYZ());
	boners.push_back(bone);
	
	if(!parent)theRoot = bone;
      }
    }
    if(theRoot){
      init(theRoot);
      useEuclideanRoot();
    }
    return (theRoot != 0);
  }

  std::vector<std::string> getObjStrings() const {
    std::vector<std::string> strings;
    
    for(int i=0; i<(int)boneByID.size(); i++){
      Bone * bone=boneByID[i];
      Mat4x4 m=bone->lmat;
      //m=bone->lanim;//Used to be lanim and then undo bone->len transform because parms were zero'd
      //m.trans_right(0,-bone->len,0);
      
      Quaternion q=Quaternion::fromMatrix(m);
      
      double px=m(0,3);
      double py=m(1,3);
      double pz=m(2,3);
      
      char line[1024];
      snprintf(line, 1024, "b %s %d %g %g %g %g %g %g %g %g",bone->name.c_str(),
	       (bone->parent)?bone->parent->id:-1,bone->len,q.v.x,q.v.y,q.v.z,q.a,px,py,pz);
      strings.push_back(std::string(line));
    }    
    return strings;
  }

  void initHumanTips();
};


//FIXME: move these guys someplace, for generating reports.
std::pair<double, Matrix> compareParamsDirect(Armature & arm, const Matrix & p1, const Matrix & p2);
std::pair<double, Matrix> compareParamsQuats(Armature & arm, const Matrix & p1, const Matrix & p2);

void compareParamsSummary(Armature & arm, PoseKeys & k1, PoseKeys & k2);


#endif
