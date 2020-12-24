#include <GL/glew.h>
#include "mesh.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <nmath/matrix.h>
#include <assert.h>
#include <nmath/vec4.h>
#include <nmath/mat4x4.h>
#include <map>

using nacb::Vec4f;
using nacb::Vec3f;
using nacb::Vec2f;
using nacb::Mat4x4;
using std::map;




void drawTC(Mesh * mesh, const nacb::Vec3d & ccd){
  glColor3f(1,1,1);
  nacb::Vec3f cc = ccd;
  glBegin(GL_TRIANGLES);
  for(int i=0; i<(int)mesh->tris.size(); i++){
    double ds[3]; // double ds...hahaha
    int nless = 0;

    for(int k=0; k<3; k++){
      Vec3f v = cc - mesh->vert[mesh->tris[i].vi[k]];
      v.normalize();
      double d = v.dot(mesh->norm[mesh->tris[i].ni[k]]);

      ds[k] = d;
      nless += (d <= 0);
    }
    
    // If all three are behind don't draw...this didn't do what I expected.
    if (nless >= 2) {
      continue;
    }

    for (int k=0; k<3; k++) {
      glColor4f(1,1,1,std::max(0.0, ds[k]));
      glMultiTexCoord3fv(GL_TEXTURE1, mesh->vert[mesh->tris[i].vi[k]].data);
      glTexCoord3fv(mesh->vert[mesh->tris[i].vi[k]].data);
      glVertex2fv(mesh->tvert[mesh->tris[i].tci[k]].data);
    }
  }
  glEnd();
}

//Used in warp texture
void depthRange(Mesh * mesh, const double * zdir, 
		double & zmin, double & zmax){
  zmin =  10e10;
  zmax = -10e10;

  for(int i=0; i<(int)mesh->vert.size(); i++){
    const Vec3f & v =  mesh->vert[i];
    double z = zdir[0]*v.x + zdir[1]*v.y + zdir[2]*v.z + zdir[3];
    if (isnan(z) || isinf(z)) continue;
    zmin = std::min(z, zmin);
    zmax = std::max(z, zmax);
  }
  zmin -= 5e-2;
  zmax += 5e-2;
}


static void normalize(Mesh::weight_vector_t & wts){
  double sum = 0;
  for(int i=0; i<(int)wts.size(); i++){
    sum += wts[i].weight;
  }
  if(sum>1e-10){
    for(int i=0; i<(int)wts.size(); i++)
      wts[i].weight *= (1.0/sum);    
  }
}

static void append(Mesh::weight_vector_t & dest, const Mesh::weight_vector_t & src){
  for(int i=0; i<(int)src.size(); i++){
    bool fnd = false;
    for(int j=0; j<(int)dest.size() && !fnd; j++){
      if(dest[j].bone == src[i].bone){
	dest[j].weight += src[i].weight;
	fnd = true;
	break;
      }
    }
    if(!fnd){
      dest.push_back(src[i]);
    }
  }
}

bool operator==(const vector<Mesh::bone_weight_t> & v1,
		const vector<Mesh::bone_weight_t> & v2){
  return &v1==&v2;
}

Mesh::Mesh(const nacb::Imagef & coords, const nacb::Image8 & mask){
  assert(coords.nchannels == 3);
  assert(mask.w == coords.w && mask.h == coords.h);

  vinds = 0;
  vinds_size = 0;
  vbo = 0;

  nacb::Image32 inds(mask.w, mask.h, 1);
  inds = -1;
  for(int y=0; y<coords.h; y++){
    for(int x=0; x<coords.w; x++){
      if(mask(x, y)){
	inds(x, y) = restVert.size();
	restVert.push_back(Vec3f(coords(x, y, 0),
				 coords(x, y, 1),
				 coords(x, y, 2)));
	tvert.push_back(Vec3f(float(x)/coords.w, float(y)/coords.h, 0.0));
      }
    }
  }

  for(int y=0; y<coords.h-1; y++){
    for(int x=0; x<coords.w-1; x++){
      if(mask(x, y) && mask(x+1, y) && mask(x+1, y+1)){
	int ix[3] = {(int)inds(x, y), (int)inds(x+1, y+1), (int)inds(x+1, y)};	
	tris.push_back(Mesh::Tri(ix, ix, ix));
      }
      if(mask(x, y) && mask(x+1, y+1) && mask(x, y+1)){
	int ix[3] = {(int)inds(x, y), (int)inds(x, y+1), (int)inds(x+1, y+1)};	
	tris.push_back(Mesh::Tri(ix, ix, ix));
      }
    }
  }

  bone_weights = bone_weights_t(restVert.size());

  vert = restVert;
  
  initNormals();
}

void Mesh::initVertexElements(){
  if(!vinds || vinds_size != (int)tris.size()*9)
  {
    delete [] vinds;
    vinds = new float[tris.size()*9];
    vinds_size = tris.size()*9;
  }
  for(int i=0; i<(int)tris.size(); i++){
    for(int k=0; k<3; k++)
      memcpy(vinds+9*i+3*k, vert[tris[i].vi[k]].data, sizeof(float)*3);
  }
  /*
  if(!vbo){
    glGenBuffers(1, &vbo);
  }
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vinds_size, vinds, GL_DYNAMIC_DRAW_ARB);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  */
}

void Mesh::drawVertexElements(){
  //assert(tris.size()*9 == vinds_size);
  glEnableClientState(GL_VERTEX_ARRAY);

  if(vbo){
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexPointer(3, GL_FLOAT, 0, 0);
  }
  else
    glVertexPointer(3, GL_FLOAT,0, vinds);
  
  glDrawArrays(GL_TRIANGLES, 0, tris.size()*3);
  
  if(vbo)
    glBindBuffer(GL_ARRAY_BUFFER, 0);

  glDisableClientState(GL_VERTEX_ARRAY);
}

void Mesh::drawWithTC(const vector<Vec3f> & tc){
  glBegin(GL_QUADS);
  for(int i=0;i<(int)poly.size();i++){
    for(int j=0;j<4;j++){
      glTexCoord3fv(tc[poly[i].vi[j]].data);
      glNormal3fv(norm[poly[i].ni[j]].data);
      glVertex3fv(vert[poly[i].vi[j]].data);
    }
  }
  glEnd();
  
  glBegin(GL_TRIANGLES);
  for(int i=0;i<(int)tris.size();i++){
    for(int j=0;j<3;j++){
      glTexCoord3fv(tc[tris[i].vi[j]].data);
      glNormal3fv(norm[tris[i].ni[j]].data);
      glVertex3fv(vert[tris[i].vi[j]].data);
    }
  }
  glEnd();
}

void Mesh::draw(bool fancy){
  glPushMatrix();
  
  glBegin(GL_QUADS);
  for(int i=0;i<(int)poly.size();i++){
    for(int j=0;j<4;j++){
      glTexCoord3fv(tvert[poly[i].tci[j]].data);
      glNormal3fv(norm[poly[i].ni[j]].data);
      glVertex3fv(vert[poly[i].vi[j]].data);
    }
  }
  glEnd();
  
  glBegin(GL_TRIANGLES);
  for(int i=0;i<(int)tris.size();i++){
    for(int j=0;j<3;j++){
      glTexCoord3fv(tvert[tris[i].tci[j]].data);
      glNormal3fv(norm[tris[i].ni[j]].data);
      glVertex3fv(vert[tris[i].vi[j]].data);
    }
  }
  glEnd();
  glPopMatrix();

  if(fancy){
    GLUquadric * quad=gluNewQuadric();
    for(int i=0;i<(int)spheres.size();i++){
      bool colliding=false;
      for(int j=0;j<(int)spheres[i].interacts.size();j++){
	int other=spheres[i].interacts[j];
	double d=(spheres[i].c-spheres[other].c).len();
	if(d<=spheres[i].r+spheres[other].r){
	  colliding=true;
	  break;
	}
      }
      if(colliding)
	glColor4f(0,0,1,0.4);
      else
	glColor4f(1,1,1,0.4);

      glPushMatrix();
      glTranslatef(spheres[i].c.x,spheres[i].c.y,spheres[i].c.z);
      gluSphere(quad,spheres[i].r,32,32);
      glPopMatrix();

      /*
	glBegin(GL_LINES);
	for(int j=0;j<(int)spheres[i].interacts.size();j++){
	glVertex3fv(spheres[i].c.data);
	glVertex3fv(spheres[spheres[i].interacts[j]].c.data);
	}
	glEnd();
      */
    }
    glColor3f(1,1,1);
    gluDeleteQuadric(quad);
  }
}

void Mesh::initNormalsAutoSmooth(double degreeThresh) {
  // Each vertex in each faces gets its own normal.
  norm = std::vector<Vec3f>(tris.size()*3, Vec3f(0, 0, 0));
  
  vector<Vec3f> faceNorm(tris.size(), Vec3f(0, 0, 0));
  vector<vector<int> > vertFaces(vert.size());
  for (int i=0; i<(int)tris.size(); i++) {
    Vec3f e1 = vert[tris[i].vi[1]] - vert[tris[i].vi[0]];
    Vec3f e2 =vert[tris[i].vi[2]] - vert[tris[i].vi[1]];
    Vec3f n = e1.cross(e2);

    n.normalize();

    faceNorm[i] = n;

    for (int k=0; k<3; k++) {
      vertFaces[tris[i].vi[k]].push_back(i);
    }
  }

  double thresh = cos(degreeThresh * M_PI / 180.0);
  for (int i=0; i<(int)tris.size(); i++) {
    for (int k=0; k<3; k++) {
      int vi = tris[i].vi[k];
      
      Vec3f n(0, 0, 0);
      for (int h=0; h<(int)vertFaces[vi].size(); h++) {
	int f = vertFaces[vi][h];
	if (faceNorm[f].dot(faceNorm[i]) >= thresh) {
	  n += faceNorm[f];
	}
      }
      n.normalize();
      norm[3*i + k] = n;
      tris[i].ni[k] = 3*i + k;
    }
  }
}

void Mesh::initNormalsFlat(){
  norm = std::vector<Vec3f>(poly.size() + tris.size());

  for(int i=0; i<(int)poly.size(); i++){
    Vec3f e1 = vert[poly[i].vi[1]]-vert[poly[i].vi[0]];
    Vec3f e2 = vert[poly[i].vi[2]]-vert[poly[i].vi[1]];
    Vec3f n = e1.cross(e2);
    n.normalize();

    for(int j=0; j<4; j++){
      poly[i].ni[j] = i;
      norm[i] = n;
    }	
  }

  for(int i=0;i<(int)tris.size();i++){
    Vec3f e1=vert[tris[i].vi[1]]-vert[tris[i].vi[0]];
    Vec3f e2=vert[tris[i].vi[2]]-vert[tris[i].vi[1]];
    Vec3f n=e1.cross(e2);

    for(int j=0; j<3; j++){
      tris[i].ni[j] = i + poly.size();
      norm[i + poly.size()] = n;
    }	
  }
}

void Mesh::initNormals(){
  norm.clear();
  for(int i=0;i<(int)vert.size();i++)
    norm.push_back(Vec3f(0,0,0));
  for(int i=0;i<(int)poly.size();i++){
    Vec3f e1=vert[poly[i].vi[1]]-vert[poly[i].vi[0]];
    Vec3f e2=vert[poly[i].vi[2]]-vert[poly[i].vi[1]];
    Vec3f n=e1.cross(e2);
    for(int j=0;j<4;j++){
      poly[i].ni[j]=poly[i].vi[j];
      norm[poly[i].ni[j]]+=n;
    }	
  }
  for(int i=0;i<(int)tris.size();i++){
    Vec3f e1=vert[tris[i].vi[1]]-vert[tris[i].vi[0]];
    Vec3f e2=vert[tris[i].vi[2]]-vert[tris[i].vi[1]];
    Vec3f n=e1.cross(e2);
    for(int j=0;j<3;j++){
      tris[i].ni[j]=tris[i].vi[j];
      norm[tris[i].ni[j]]+=n;
    }	
  }
  for(int i=0;i<(int)norm.size();i++)
    norm[i].normalize();
}

int Mesh::saveobj(const char * fname,bool append, bool deformed){
  FILE * file=fopen(fname,append?"a":"w");
  if(!file)return 0;

  if(deformed){
    for(int i=0;i<(int)vert.size();i++)
      fprintf(file,"v %g %g %g\n",vert[i].x,vert[i].y,vert[i].z);
  }
  else {
    for(int i=0;i<(int)restVert.size();i++)
      fprintf(file,"v %g %g %g\n",restVert[i].x,restVert[i].y,restVert[i].z);
  }
  for(int i=0;i<(int)tvert.size();i++){
    fprintf(file,"vt %g %g\n",tvert[i].x,tvert[i].y);
  }
  for(int i=0;i<(int)tris.size();i++){
    if(tvert.size())
      fprintf(file,"f %d/%d %d/%d %d/%d\n",
	      tris[i].vi[0]+1,
	      tris[i].tci[0]+1,
	      tris[i].vi[1]+1,
	      tris[i].tci[1]+1,
	      tris[i].vi[2]+1,
	      tris[i].tci[2]+1);
    else
      fprintf(file,"f %d %d %d\n",tris[i].vi[0]+1,
	      tris[i].vi[1]+1,
	      tris[i].vi[2]+1);
  }

  if(bone_weights.size()){
    for(int i=0;i<(int)restVert.size();i++){
      for(int k=0;k<(int)bone_weights[i].size();k++){
	fprintf(file,"vw %d %d %g\n",i+1,bone_weights[i][k].bone,
		bone_weights[i][k].weight);
      }
    }
  }

  for(int i=0;i<(int)spheres.size();i++){
    fprintf(file,"sp %g %g %g %g\n",
	    spheres[i].c.x,spheres[i].c.y,spheres[i].c.z,
	    spheres[i].r);
  }
  for(int i=0;i<(int)spheres.size();i++){
    for(int k=0;k<(int)spheres[i].weights.size();k++){
      fprintf(file,"sw %d %d %g\n",i+1,spheres[i].weights[k].bone,
	      spheres[i].weights[k].weight);
    }
  }
  for(int i=0;i<(int)spheres.size();i++){
    for(int k=0;k<(int)spheres[i].interacts.size();k++){
      int other=spheres[i].interacts[k];
      //only output once, they are symmetric
      if(i<other){
	fprintf(file,"si %d %d\n",i+1,other+1);
      }
    }
  }

  fclose(file);
  return 1;
}

int Mesh::loadobj(const char * fname,
		  vector<std::string> * unused){
  vert.clear();
  tvert.clear();
  norm.clear();
  poly.clear();
  tris.clear();
  spheres.clear();
  bone_weights.clear();
  vector<int> boneParents;

  FILE * file=fopen(fname,"r");
  if(!file)return 0;
  char line[1024];
  int nonorm=0;
  bool got_bone_weights=0;

  while(fgets(line,1024,file)!=0){
#ifdef USE_BONES
    if(line[0]=='b'){
      int parenti;
      double len;
      double qx,qy,qz,qw;
      double px,py,pz;
      char name[1024];
      if(10==sscanf(line+2,"%s %d %lf %lf %lf %lf %lf %lf %lf %lf",
		    name,
		    &parenti,&len,&qx,&qy,&qz,&qw,&px,&py,&pz)){
	printf("got bone \"%s\" %d %f %f %f  %f %f %f %f\n",name,parenti,qx,qy,qz,qw,px,py,pz);
	quat q;
	q[0]=qw;
	q[1]=qx;
	q[2]=qy;
	q[3]=qz;
	bones.push_back(new Bone(0,name,Vec3f(px,py,pz),q,len));
	bones[bones.size()-1]->setIndex(bones.size()-1);
	boneParents.push_back(parenti);
      }
    }
#endif
    if((line[0]=='f'||line[0]=='F') && line[1]==' '){
      int v[4]={1,1,1,1};
      int t[4]={1,1,1,1};
      int n[4]={1,1,1,1};
      if(12==sscanf(line+2,"%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",
		    v,t,n,v+1,t+1,n+1,v+2,t+2,n+2,v+3,t+3,n+3)){
	for(int i=0;i<4;i++)v[i]--,t[i]--,n[i]--;

	poly.push_back(Face(v,n,t));
      }
      else if(8==sscanf(line+2,"%d/%d %d/%d %d/%d %d/%d",
			v,t,v+1,t+1,v+2,t+2,v+3,t+3)){
	for(int i=0;i<4;i++)v[i]--,t[i]--,n[i]--;
	poly.push_back(Face(v,v,t));
	nonorm=1;
      }
      else if(4==sscanf(line+2,"%d %d %d %d",
			v,v+1,v+2,v+3)){
	for(int i=0;i<4;i++)v[i]--,t[i]--,n[i]--;
	poly.push_back(Face(v,v,t));
	nonorm=1;
      }
      else{
	if(9==sscanf(line+2,"%d/%d/%d %d/%d/%d %d/%d/%d",
		     v,t,n,v+1,t+1,n+1,v+2,t+2,n+2)){
	  for(int i=0;i<3;i++)v[i]--,t[i]--,n[i]--;
	  tris.push_back(Tri(v,n,t));
	}
	else if(6==sscanf(line+2,"%d/%d %d/%d %d/%d",
			  v,t,v+1,t+1,v+2,t+2)){
	  for(int i=0;i<3;i++)v[i]--,t[i]--,n[i]--;
	  tris.push_back(Tri(v,v,t));
	  nonorm=1;
	}
	else if(3==sscanf(line+2,"%d %d %d",
			  v,v+1,v+2)){
	  for(int i=0;i<3;i++)v[i]--,t[i]--,n[i]--;
	  tris.push_back(Tri(v,v,t));
	  nonorm=1;
	}
      }
    }
    else if(line[0]=='v' && line[1]=='t'){
      Vec3f tc(0,0,1);
      if(2==sscanf(line+3,"%f %f",tc.data,tc.data+1)){
	tvert.push_back(tc);
      }
      else{
	printf("error reading tvert\n");
	return 0;
      }
    }
    else if(line[0]=='v' && line[1]=='n'){
      Vec3f vec;
      if(3==sscanf(line+3,"%f %f %f",&vec.x,&vec.y,&vec.z)){
	vec.normalize();
	norm.push_back(vec);
      }
      else{
	printf("error reading vn\n");
	return 0;
      }
    }
    else if(line[0]=='v' && line[1]==' '){
      Vec3f vec;
      if(3==sscanf(line+2,"%f %f %f",&vec.x,&vec.y,&vec.z)){
	vert.push_back(vec);
	bone_weights.push_back(vector<bone_weight_t>());
      }
      else{
	printf("error reading vn\n");
	return 0;
      }
    }
    else if(tolower(line[0])=='v' && tolower(line[1])=='w'){
      int vert;
      int bone;
      float weight;
      if(3==sscanf(line+3,"%d %d %f",&vert,&bone,&weight)){
	//printf("got bone weight %d,%d %f\n",vert,bone,weight);
	vert=vert-1;//just to be like .obj
	bone_weights[vert].push_back(bone_weight_t(bone,weight));
	got_bone_weights=true;
      }
    }
    else if(tolower(line[0])=='s' && tolower(line[1])=='p'){
      Vec3f p;
      float r;
      if(4==sscanf(line+3,"%f %f %f %f",&(p.x),&(p.y),&(p.z),&r)){
	addSphere(p,r,vector<bone_weight_t>());
      }
    }
    else if(tolower(line[0])=='s' && tolower(line[1])=='w'){
      int spherei;
      int bone;
      float weight;
      if(3==sscanf(line+3,"%d %d %f",&spherei,&bone,&weight)){
	//printf("got bone weight %d,%d %f\n",vert,bone,weight);
	spherei=spherei-1;//just to be like .obj
	if(spherei>=(int)spheres.size()){
	  fprintf(stderr,"ERROR: sphere weights are going to be wrong, fix this loader\n");
	}
	else
	  spheres[spherei].weights.push_back(bone_weight_t(bone,weight));
      }
    }
    else if(tolower(line[0])=='s' && tolower(line[1])=='i'){
      int sphere1;
      int sphere2;
      if(2==sscanf(line+3,"%d %d",&sphere1,&sphere2)){
	sphere1--;
	sphere2--;
	if(sphere1<0 || sphere2<0  ||
	   sphere1>=(int)spheres.size() ||
	   sphere2>=(int)spheres.size()){
	  fprintf(stderr,"ERROR: sphere interaction out of bounds, must come after sphere definitions");
	}
	else{
	  printf("got sphere interaction %d %d\n",sphere1,sphere2);
	  spheres[sphere1].interacts.push_back(sphere2);
	  spheres[sphere2].interacts.push_back(sphere1);
	}
      }
    }
    else if(unused){
      unused->push_back((std::string)line);
    }
  }
  fclose(file);

  if(!tvert.size()){
    tvert.push_back(Vec3f(0,0,1));
    for(int i=0;i<(int)poly.size();i++){
      for(int j=0;j<4;j++){
	if(poly[i].tci[j]>=(int)tvert.size()){
	  printf("texture coordinate index out of bounds:%d \n",poly[i].tci[j]);
	}
	if(poly[i].tci[j]<0)poly[i].tci[j]=0;
      }
    }
    for(int i=0;i<(int)tris.size();i++){
      for(int j=0;j<3;j++){
	if(tris[i].tci[j]>=(int)tvert.size()){
	  printf("texture coordinate index out of bounds:%d \n",tris[i].tci[j]);
	}
	if(tris[i].tci[j]<0)tris[i].tci[j]=0;
      }
    }
  }
  restVert=vert;

  if(got_bone_weights){
    for(int i=0;i<(int)bone_weights.size();i++){
      float sum=0;
      for(int j=0;j<(int)bone_weights[i].size();j++){
	sum+=bone_weights[i][j].weight;
      }
      if(sum<=1e-6)printf("no weight for vert %d\n",i);
      for(int j=0;j<(int)bone_weights[i].size();j++){
	bone_weights[i][j].weight/=sum;
      }
    }
  }
  else bone_weights.clear();

  for(int i=0;i<(int)spheres.size();i++){
    spheres[i].cleanInteracts();
  }

  if(nonorm){
    //printf("there were no normals\n");
    initNormals();
  }
  initNormals();
  return 1;
}

template <class T>
void Mesh::animate(T tforms){
  //This would cause a crashie
  if(bone_weights.size()<restVert.size())return;

  for(int i=0; i<(int)restVert.size(); i++){
    const Vec3f & v=restVert[i];
     Vec3f sv(0,0,0);
     double sum=0;
    for(int j=0;j<(int)bone_weights[i].size();j++){
      sv+=(tforms[bone_weights[i][j].bone].transformPoint(v))*bone_weights[i][j].weight;
       sum+=bone_weights[i][j].weight;
    }
    if(fabs(sum-1.0)>1e-6)sv*=(1.0/sum);
    vert[i]=sv;
  }
  initNormals();

  for(int i=0;i<(int)spheres.size();i++){
    const Vec3f & v=spheres[i].restc;
    Vec3f sv(0,0,0);
    double sum=0;
    for(int j=0;j<(int)spheres[i].weights.size();j++){
      sv+=(tforms[spheres[i].weights[j].bone].transformPoint(v))*spheres[i].weights[j].weight;
      sum+=spheres[i].weights[j].weight;
    }
    if(fabs(sum-1.0)>1e-6)sv*=(1.0/sum);
    spheres[i].c=sv;
  }
}

Mesh Mesh::subdivide() const{
  if(poly.size())
    fprintf(stderr, "WARNING: subdivide does not handle polys yet\n");

  map<int, int> * adj = new map<int,int>[vert.size()];
  map<int, int> * tadj = new map<int,int>[tvert.size()];
  vector<Vec3f> newVert;
  vector<Tri> newTris;
  bone_weights_t newWeights;
  vector<Vec3f> newTC;
  bool has_tc = tvert.size();

  newVert = restVert;
  newWeights = bone_weights;
  newTC = tvert;

  for(int i=0; i<(int)tris.size(); i++){
    Vec3f avg(0,0,0);
    Vec3f tavg(0,0);
    int vcent = newVert.size();
    int tvcent = newTC.size();

    newVert.push_back(avg);
    newWeights.push_back(weight_vector_t());

    if(has_tc){
      newTC.push_back(tavg);
    }
    for(int k=0; k<3; k++){
      int  v1 = tris[i].vi[k];
      int  v2 = tris[i].vi[(k+1)%3];
      int tv1 = tris[i].tci[k];
      int tv2 = tris[i].tci[(k+1)%3];
     
      avg += restVert[v1];

      if(!adj[v1].count(v2)){
	adj[v1][v2] = newVert.size();
	adj[v2][v1] = newVert.size();
	newVert.push_back((restVert[v1]+restVert[v2])*0.5);

	weight_vector_t wts = bone_weights[v1];
	append(wts, bone_weights[v2]);
	normalize(wts);
	newWeights.push_back(wts);
      }

      if(has_tc)
	tavg += tvert[tv1];
      
      if(has_tc && !tadj[tv1][tv2]){
	tadj[tv1][tv2] = newTC.size();
	tadj[tv2][tv1] = newTC.size();

	newTC.push_back((tvert[tv1]+tvert[tv2])*0.5);
      }

      int vi[3];
      int tci[3] = {0,0,0};
      vi[0] = v1;
      vi[1] = adj[v1][v2];
      vi[2] = vcent;

      if(has_tc){
	tci[0] = tv1;
	tci[1] = tadj[tv1][tv2];
	tci[2] = tvcent;
      }
      newTris.push_back(Tri(vi, vi, tci));

      vi[0] = adj[v1][v2];
      vi[1] = v2;
      vi[2] = vcent;

      if(has_tc){
	tci[0] = tadj[tv1][tv2];
	tci[1] = tv2;
	tci[2] = tvcent;
      }
      newTris.push_back(Tri(vi, vi, tci));

      append(newWeights[vcent], bone_weights[v1]);
    }
    if(has_tc)
      newTC[tvcent] = tavg*(1.0/3.0);
    normalize(newWeights[vcent]);
    newVert[vcent] = avg*(1.0/3.0);
  }
  delete [] adj;

  Mesh newMesh;
  newMesh.restVert = newVert;
  newMesh.vert = newVert;
  newMesh.tris = newTris;
  newMesh.bone_weights = newWeights;
  newMesh.initNormals();
  newMesh.spheres = spheres;
  newMesh.tvert = newTC;

  return newMesh;
}


Mesh Mesh::keepVert(const std::vector<int> & vi) const {
  std::vector<int> vmap(vert.size(), -1);
  std::vector<int> tcmap(tvert.size(), -1);
  std::vector<Vec3f>  vnew;
  std::vector<Vec3f> tcnew;
  bone_weights_t      wnew;
  std::vector<Tri>  trinew;

  for(int i=0; i<(int)vi.size(); i++){
    vmap[vi[i]] = i;
    vnew.push_back(vert[vi[i]]);
    wnew.push_back(bone_weights[vi[i]]);
  }
  
  for(int i=0; i<(int)tris.size(); i++){
    bool good = true;
    for(int j=0; j<3; j++){
      good &= (vmap[tris[i].vi[j]]>=0);
    }
    if(good){
      int vi[3], ti[3], ni[3]={0,0,0};
      for(int j=0; j<3; j++){
	vi[j] = vmap[tris[i].vi[j]];
	if(tcmap[tris[i].tci[j]]<0){
	  tcmap[tris[i].tci[j]] = tcnew.size();
	  tcnew.push_back(tvert[tris[i].tci[j]]);
	}
	ti[j] = tcmap[tris[i].tci[j]];
      }
      trinew.push_back(Tri(vi, ni, ti));
    }
  }
  //FIXME: doesn't copy spheres.
  Mesh m;
  m.vert = vnew;
  m.restVert = vnew;
  m.bone_weights = wnew;
  m.tvert = tcnew;
  m.tris = trinew;
  m.initNormals();
  return m;
}


void Mesh::fillSimpleHoles(){
  std::map<int, int> edges;
  for(int i=0; i<(int)tris.size(); i++){
    for(int j=0; j<3; j++){
      int   e = tris[i].vi[j] * vert.size() + tris[i].vi[(j+1)%3];
      edges[e] = tris[i].vi[j] + tris[i].vi[(j+1)%3] * vert.size();
    }
  }
  std::map<int, int>::iterator f = edges.begin(), end = edges.end();
  int * missing = new int[vert.size()];
  int * visited = new int[vert.size()];
  int nmiss = 0;
  memset(missing, 0xFF, sizeof(int)*vert.size());
  memset(visited, 0xFF, sizeof(int)*vert.size());

  while(f!=end){
    if(!edges.count(f->second)){
      int v1 = f->first/vert.size();
      int v2 = f->first % vert.size();
      
      missing[v2] = v1;
      nmiss++;
    }
    f++;
  }
  if(nmiss){
    for(int i=0; i<(int)vert.size(); i++){
      if(missing[i] == -1)continue;
      
      int finger = i;
      int cnt = 0;
      
      while(finger>=0 && visited[finger] != i){
	visited[finger] = i;
	finger = missing[finger];
	cnt++;
      }

      if(finger == i && cnt){
	printf("found a hole of size: %d, filling it\n", cnt);

	std::vector<int> vfill;

	finger = i;
	while(finger>=0){
	  vfill.push_back(finger);

	  missing[finger] = -1;
	  finger = missing[finger];
	  cnt++;
	}

	for(int i=0; i<int(vfill.size())-2; i++){
	  int vi[3] = {vfill[0], vfill[i+1], vfill[i+2]};
	  int ni[3] = {0, 0, 0};
	  tris.push_back(Tri(vi, ni));
	}
      }
    }
  }

  delete [] visited;
  delete [] missing;
  initNormals();
  return;
}


void intersectRayWithMesh(const Mesh & mesh,
			  const Vec3f & o,const Vec3f & k,
			  double & mindist,
			  double & maxdist,
			  int & triangleIndex,
			  Vec3f & bary, 
			  bool returnNear){
  mindist=10e10;
  maxdist=-10e10;
  for(int i=0;i<(int)mesh.tris.size();i++){
    Vec3f v0=mesh.restVert[mesh.tris[i].vi[0]];
    Vec3f v1=mesh.restVert[mesh.tris[i].vi[1]];
    Vec3f v2=mesh.restVert[mesh.tris[i].vi[2]];
    Vec3f e1=v1-v0;
    Vec3f e2=v2-v0;
    Vec3f n=e1.cross(e2);
    n.normalize();
    double d=-v0.dot(n);
    double t=(-d-n.dot(o))/(n.dot(k));
    Vec3f isect=o+k*t-v0;
    
    double A=e1.dot(e1);
    double B=e2.dot(e1);
    double C=B;
    double D=e2.dot(e2);
    
    double det=A*D-B*C;
    nacb::Vec2d ixy(e1.dot(isect), e2.dot(isect));
    double a0=( D*ixy.x-B*ixy.y)/det;
    double a1=(-C*ixy.x+A*ixy.y)/det;

    //printf("should be zero:%f\n",(e1*a0+e2*a1+v0-(cc+k*t)).len());
    if(a0>=0 && a1>=0 && a0+a1<=1){
      if((returnNear  && t<mindist) ||
	 (!returnNear && t>maxdist)){
	triangleIndex = i;
	bary = Vec3f(1.0-(a0+a1), a0, a1);
      }
      mindist=std::min(mindist,t);
      maxdist=std::max(maxdist,t);
    }
    if(a0>=-1e-5 && a1>=-1e-5 && (a0 + a1 <= 1.0 + 2e-5)){
      printf("close: %d\n", i);
    }
  }
}


void intersectRayWithMesh(const Mesh & mesh,const Vec3f & o,const Vec3f & k,
			  double & mindist,
			  double & maxdist){

  int  index;
  Vec3f bary;
  intersectRayWithMesh(mesh, o, k, mindist, maxdist, index, bary, true);
}



/**
   Return the bary-centric coordinate image for the mesh 
   of size tw, th (zeros where there is no triangle).

   If timage!=null then fill it with the triangle indices (-1 
   where there is no triangle).
*/
nacb::Imagef getBaryImage(Mesh & mesh, int tw, int th, 
			  bool useTexCoords,
			  nacb::Image32 * timage){
  nacb::Imagef bary(tw, th, 4);
  nacb::Image8 trgb(tw, th, 3);

  if(timage != 0){
    *timage = nacb::Image32(tw, th, 1);
  }

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);

  glViewport(0, 0, tw, th);

  //Set-up texture matrices
  if(useTexCoords){
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0,1,0,1,-1,1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
  }
  glDisable(GL_CULL_FACE);

  glBegin(GL_TRIANGLES);
  for(int i=0; i<(int)mesh.tris.size(); i++){
    for(int k=0; k<3; k++){
      glColor3f(k==0, k==1, k==2);
      if(useTexCoords)glVertex2fv(mesh.tvert[mesh.tris[i].tci[k]].data);
      else glVertex3fv(mesh.vert[mesh.tris[i].vi[k]].data);
    }
  }
  glEnd();
  
  glReadPixels(0, 0, tw, th, bary.channelToGL(), bary.typeToGL(), bary.data);

  if(timage){
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glBegin(GL_TRIANGLES);
    for(int i=0; i<(int)mesh.tris.size(); i++){
      for(int k=0; k<3; k++){
	float r = double(i % 256)/255.0;
	float g = double((i >> 8) % 256)/255.0;
	float b = double((i >> 16) % 256)/255.0;
	glColor3f(r, g, b);

	if(useTexCoords)glVertex2fv(mesh.tvert[mesh.tris[i].tci[k]].data);
	else glVertex3fv(mesh.vert[mesh.tris[i].vi[k]].data);
      }
    }
    glEnd();

    glReadPixels(0, 0, tw, th, trgb.channelToGL(), trgb.typeToGL(), trgb.data);
  }

  if(useTexCoords){
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
  }

  for(int y=0; y<bary.h; y++){
    for(int x=0; x<bary.w; x++){
      nacb::Vec3f co(bary(x,y,0), bary(x,y,1), bary(x,y,2));
      double len = co.x + co.y + co.z;
      if(len<1e-2){
	bary(x,y,0) = bary(x,y,1) = bary(x,y,2) = bary(x,y,3) = 0;
	if(timage)(*timage)(x,y) = uint32_t(-1);
      }
      else {
	co *= (1.0/len);

	if(timage){
	  int ti = int(trgb(x,y,0)) + (int(trgb(x,y,1))<<8) + (int(trgb(x,y,2))<<16);
	  (*timage)(x,y) = ti;
	}
	for(int k=0; k<3; k++)
	  bary(x,y,k) = co.data[k];
	bary(x,y,3) = 1.0;
      }
    }
  }
  return bary;
}


template void Mesh::animate(nacb::Mat4x4 *);
template void Mesh::animate(std::vector<nacb::Mat4x4, std::allocator<nacb::Mat4x4> >);
