#include <nmath/matrix.h>
#include <nimage/image.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include <nmisc/timer.h>
#include <vector>
#include <GL/gl.h>
#include <sys/stat.h>
#include <assert.h>
#include "utils.h"
#include "../recon_globals.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <zlib.h>

using namespace std;
using namespace nacb;
namespace bfs=boost::filesystem;


void uniform3x3(GLint loc, const Matrix & m){
  static float data[9];
  for(int i=0; i<9; i++)
    data[i] = m.data[i];
  glUniformMatrix3fv(loc, 1, true, data);
}


void applyHomographyGL(GLuint prog,
		       const Matrix & K,
		       const Matrix & Kinv,
		       const Matrix & dd,
		       const Matrix & H,
		       int w, int h,
		       double minx,
		       double maxx,
		       double miny,
		       double maxy){
  GLuint loc;
  glUseProgram(prog);

  loc = glGetUniformLocation(prog, "k");
  assert(loc>=0);
  float d[5];
  for(int i=0; i<5; i++)d[i] = dd[i];
  glUniform1fv(loc, 5, d);

  loc = glGetUniformLocation(prog, "Kinv");assert(loc>=0);
  uniform3x3(loc, K.inverse());
  
  loc = glGetUniformLocation(prog, "K");assert(loc>=0);
  uniform3x3(loc, K);
  
  loc = glGetUniformLocation(prog, "imsize");assert(loc>=0);
  glUniform2f(loc, w, h);

  ::applyHomographyGL(H, minx, maxx, miny, maxy);

  glUseProgram(0);
}

void applyHomographyGL(const Matrix & H,
		       double minx,
		       double maxx,
		       double miny,
		       double maxy){
  glPushAttrib(GL_VIEWPORT_BIT);

  Matrix Hinv = H.inverse();
  glViewport(0, 0, (int)ceil(maxx-minx), (int)ceil(maxy-miny));

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(minx, maxx, miny, maxy, -100.0, 100.0);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  glLoadIdentity();

  Matrix Huse = Matrix::eye(4,4);
  int maps[4]={0,1,-1,2};

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      if(maps[i]>=0 && maps[j]>=0)
	Huse(i,j) = Hinv(maps[i],maps[j]);
    }
  }
  Huse = Huse.transpose();
  glMultMatrixd(Huse.data);
  
  ///Huse.printMatlab("Huse");
  ///printf("range %f,%f to %f,%f\n", minx, miny, maxx, maxy);

  glBegin(GL_QUADS);
  glTexCoord4f(minx, miny, 0, 1.0);
  glVertex2f(minx, miny);

  glTexCoord4f(maxx, miny, 0, 1.0);
  glVertex2f(maxx, miny);
  
  glTexCoord4f(maxx, maxy, 0, 1.0);
  glVertex2f(maxx, maxy);

  glTexCoord4f(minx, maxy, 0, 1.0);
  glVertex2f(minx, maxy);
  glEnd();

  glMatrixMode(GL_TEXTURE);
  glPopMatrix();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glPopAttrib();
}

template <class T>
 Image<T> applyHomography(Image<T> im,
			  Matrix & H,
			  double minx,
			  double maxx,
			  double miny,
			  double maxy)
{
  Matrix Hinv=H.inverse();

  int width=(int)(maxx-minx); //ceil?
  int height=(int)(maxy-miny);
  Image<T> res(width,height,im.nchannels);
  for(int y=0;y<height;y++){
    for(int x=0;x<width;x++){
      Matrix pt(3,1);
      pt[0]=x+minx;
      pt[1]=y+miny;
      pt[2]=1;
      pt=Hinv*pt;
      pt[0]/=pt[2];
      pt[1]/=pt[2];
      int x0=(int)pt[0];
      int y0=(int)pt[1];
      double a=pt[0]-x0;
      double b=pt[1]-y0;
      if(x0>=0 && x0<im.w &&
	 y0>=0 && y0<im.h){
	int x1=x0+1;
	int y1=y0+1;
	if(x1==im.w)x1=x0;
	if(y1==im.h)y1=y0;
	for(int k=0;k<im.nchannels;k++){
	  res(x,y,k)=(T)((1.0-a)*(1.0-b)*((double)im(x0,y0,k))+
			 (a)*(1.0-b)*((double)im(x1,y0,k))+
			 (a)*(b)*((double)im(x1,y1,k))+
			 (1.0-a)*(b)*((double)im(x0,y1,k)));
	}
      }
      else{
	for(int k=0;k<im.nchannels;k++){
	  res(x,y,k)=(T)0.0;
	}
      }
    }
  }  
  return res;
}

void findBounds(Image8 & im1,Matrix & H1,
		double & minx,double & maxx)
{
  Matrix pts(3,4);
  pts.setAll(1);
  pts(0,0)=0;pts(1,0)=0;
  pts(0,1)=im1.w;pts(1,1)=0;
  pts(0,2)=im1.w;pts(1,2)=im1.h;
  pts(0,3)=0;pts(1,3)=im1.h; 

  pts=H1*pts;
  for(int j=0;j<pts.n;j++){
    pts(0,j)/=pts(2,j);
    pts(1,j)/=pts(2,j);
  }

  minx=pts.minInRow(0);
  maxx=pts.maxInRow(0);
}

void findBounds(Image8 & im1,Matrix & H1,
		Image8 & im2,Matrix & H2,
		double & minx,double & maxx,
		double & miny,double & maxy)
{
  Matrix pts(3,4);
  pts.setAll(1);
  pts(0,0)=0;pts(1,0)=0;
  pts(0,1)=im1.w;pts(1,1)=0;
  pts(0,2)=im1.w;pts(1,2)=im1.h;
  pts(0,3)=0;pts(1,3)=im1.h; 

  pts=H1*pts;
  for(int j=0;j<pts.n;j++){
    pts(0,j)/=pts(2,j);
    pts(1,j)/=pts(2,j);
  }

  minx=pts.minInRow(0);
  miny=pts.minInRow(1);
  maxx=pts.maxInRow(0);
  maxy=pts.maxInRow(1);

  pts.setAll(1);
  pts(0,0)=0;pts(1,0)=0;
  pts(0,1)=im2.w;pts(1,1)=0;
  pts(0,2)=im2.w;pts(1,2)=im2.h;
  pts(0,3)=0;pts(1,3)=im2.h; 

  pts=H2*pts;
  for(int j=0;j<pts.n;j++){
    pts(0,j)/=pts(2,j);
    pts(1,j)/=pts(2,j);
  }
  minx=std::min(pts.minInRow(0),minx);
  miny=std::min(pts.minInRow(1),miny);
  maxx=std::max(pts.maxInRow(0),maxx);
  maxy=std::max(pts.maxInRow(1),maxy);
}


nacb::Vec2d getDisparityRange(const nacb::Matrix & Pn1,
			      const nacb::Matrix & Pn2,
			      double minx1, double minx2,
			      double zmin, double zmax)
{
  Matrix K1, E1, K2, E2;
  factorProjectionMatrix(Pn1, K1, E1);
  factorProjectionMatrix(Pn2, K2, E2);

  Matrix cc2 = (E2.inverse()).getColumn(3);
  Matrix diff = E1*cc2;
  double baseline = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
  if(diff[0]<0)
    baseline *= -1;
  
  double f = (K1(0,0) + K2(0,0))/2.0;

  double d1 = -f*baseline/zmin - (minx2-minx1);
  double d2 = -f*baseline/zmax - (minx2-minx1);

  ///printf("%f %f before %f, %f\n", d1, d2, -f*baseline/zmin, -f*baseline/zmax);

  return Vec2d(std::min(d1, d2), std::max(d1, d2));
}


namespace balg=boost::algorithm;


// Faster loading a cached binary instead of a clunky .obj
// writeobj_raw and readobj_raw accomplish this.
int writeobj_raw(const char * fname, 
		  std::vector<Vec3<int> > & tris, 
		  std::vector<Vec3f> & vert){
  gzFile  file = gzopen(fname, "wb");
  int sizes[2] = {(int)vert.size(), (int)tris.size()};
  gzwrite(file, sizes, sizeof(int)*2);

  for(int i = 0; i < (int)vert.size(); i++)
    gzwrite(file, vert[i].data, sizeof(Vec3f));

  for(int i = 0; i < (int)tris.size(); i++)
    gzwrite(file, tris[i].data, sizeof(Vec3<int>));
  
  gzclose(file);
  return 1;
}


int readobj_raw(const char * fname, 
		std::vector<Vec3<int> > & tris, 
		std::vector<Vec3f> & vert){
  gzFile file = gzopen(fname, "rb");
  int sizes[2];
  gzread(file, sizes, sizeof(int)*2);

  vert = std::vector<Vec3f>(sizes[0]);
  tris = std::vector<Vec3<int> >(sizes[1]);
  
  for(int i = 0; i<(int)vert.size(); i++)
    gzread(file, vert[i].data, sizeof(Vec3f));
  
  for(int i = 0; i<(int)tris.size(); i++)
    gzread(file, tris[i].data, sizeof(Vec3<int>));
  
  gzclose(file);
  return 1;
}


void writeobj_raw(const char * fname, 
		  std::vector<Vec3<int> > & tris, 
		  std::vector<Vec3f> & vert, 
		  const Vec3f & sc, const Vec3f & tr){
  FILE * file = fopen(fname, "wb");
  int sizes[2] = {(int)vert.size(), (int)tris.size()};
  fwrite(sizes, sizeof(int), 2, file);

  for(int i=0; i<(int)vert.size(); i++)
    fwrite(vert[i].data, sizeof(Vec3f), 1, file);

  for(int i=0; i<(int)tris.size(); i++)
    fwrite(tris[i].data, sizeof(Vec3<int>), 1, file);
  
  fclose(file);
}


int readobj_and_cache(const char * fname, 
		       std::vector<Vec3<int> > & tris, 
		       std::vector<Vec3f> & vert, 
		       const Vec3f & sc, const Vec3f & tr){
  
  if(bfs::exists(std::string(fname) + ".raw")){
    if(bfs::last_write_time(std::string(fname) + ".raw") >
       bfs::last_write_time(std::string(fname))) {
      return readobj_raw((std::string(fname) + ".raw").c_str(), tris, vert);
    }
  }
  int succ = readobj(fname, tris, vert, sc, tr);

  // write out the binary cache.
  writeobj_raw((std::string(fname) + ".raw").c_str(), tris, vert);

  return succ;
}


int readobj(const char * fname, 
	    vector<Vec3<int> > & tris,
	    vector<Vec3f> & vert, 
	    const Vec3f & sc, const Vec3f & tr)
{
  nacb::StopWatch timer;

    FILE * file=fopen(fname,"r");
    char line[1024];

    if(!file){
      printf("error opening geometry.obj\n");
      return 0;
    }
    
    while(fgets(line,1024,file)){
      if(line[0]=='f'||
	 line[0]=='F'){
	Vec3<int> inds;
	if(3 != sscanf(line+2,"%d %d %d",inds.data,inds.data+1,inds.data+2)){
	  int d1,d2,d3;
	  sscanf(line+2, "%d/%d %d/%d %d/%d", inds.data, &d1, inds.data+1, &d2, inds.data+2, &d3);
	}
	inds.x--;
	inds.y--;
	inds.z--;
	tris.push_back(inds);
      }
      else if(line[1]==' '&&
	      (line[0]=='v'||line[0]=='V')){
	Vec3f v;
	sscanf(line+2,"%f %f %f",v.data,v.data+1,v.data+2);
		
	/*
	std::vector<std::string> splits;
	std::string src = std::string(line + 2);
	balg::split(splits, src, balg::is_any_of(" "));

	if(splits.size() == 3) {
	  v.x = atof(splits[0].c_str());
	  v.y = atof(splits[1].c_str());
	  v.z = atof(splits[2].c_str());
	}
	*/
	vert.push_back(Vec3f(v.x*sc.x+tr.x,v.y*sc.y+tr.y,v.z*sc.z+tr.z));
      }
    }
    fclose(file);

    printf("Reading took: %lf\n", double(timer));

    return 1;
}


int writeobj(const char * fname, 
	     std::vector<Vec3<int> > & tris, 
	     std::vector<nacb::Vec3f> & vert){
  FILE * file = fopen(fname, "w");
  if(!file)return 0;
  
  for(int i=0; i<(int)vert.size(); i++){
    fprintf(file, "v %f %f %f\n", vert[i].x, vert[i].y, vert[i].z);
  }
  
  for(int i=0; i<(int)tris.size(); i++){
    fprintf(file, "f %d %d %d\n", tris[i].x+1, tris[i].y+1, tris[i].z+1);
  }
  
  fclose(file);
  return 1;
}


/*
  From:
  A compact algorithm for rectification of stereo pairs.
  Fusiello et al. 2000.
  function art in the matlab code
*/
void rectify(Matrix & intr1, Matrix & extr1,
	    Matrix & intr2, Matrix & extr2,
	    Matrix & T1,Matrix & T2,Matrix & Pn1,Matrix & Pn2){
  int first3[]={0,1,2};
  Matrix rot1(3,3);
  Matrix rot2(3,3);
  Matrix t1(3,1),t2(3,1);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      rot1(i,j)=extr1(i,j);
      rot2(i,j)=extr2(i,j);
    }
    t1[i]=extr1(i,3);
    t2[i]=extr2(i,3);
  }
  Matrix proj = Matrix::eye(3,4);
  Matrix Po1=intr1*proj*extr1;
  Matrix Po2=intr2*proj*extr2;
  //Matrix c1= (Po1.getCols(first3,3).inverse()*Po1.getColumn(3))*-1;
  //Matrix c2= (Po2.getCols(first3,3).inverse()*Po2.getColumn(3))*-1;

  Matrix c1=(rot1.transpose())*t1*-1;
  Matrix c2=(rot2.transpose())*t2*-1;

  Vec3d oldx(rot1(0,0),rot1(0,1),rot1(0,2));
  Vec3d oldy(rot1(1,0),rot1(1,1),rot1(1,2));
  Vec3d oldz(rot1(2,0),rot1(2,1),rot1(2,2));
  Vec3d oldz_avg = oldz + Vec3d(rot2(2,0),rot2(2,1),rot2(2,2));
  oldz_avg.normalize();

  Vec3d v1(c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2]);
  //Make sure it points in the same direction..
  if(oldx.dot(v1)<0)v1*=-1;
  Vec3d v2=(oldz_avg).cross(v1);
  //Make sure y points in same direction
  if(oldy.dot(v2)<0)v2*=-1;
  Vec3d v3=v1.cross(v2);
  Matrix R(3,3);
  v1.normalize();
  v2.normalize();
  v3.normalize();
  R(0,0)=v1.x;R(0,1)=v1.y;R(0,2)=v1.z;
  R(1,0)=v2.x;R(1,1)=v2.y;R(1,2)=v2.z;
  R(2,0)=v3.x;R(2,1)=v3.y;R(2,2)=v3.z;
  Matrix newintr=(intr1+intr2);
  newintr*=0.5;
  newintr(0,1)=0;//remove skew

  Matrix t1n=(R*c1)*-1;
  Matrix t2n=(R*c2)*-1;

  Matrix extr1n=Matrix(3,4);
  Matrix extr2n=Matrix(3,4);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      extr1n(i,j)=R(i,j);
      extr2n(i,j)=R(i,j);
    }
  }
  for(int i=0;i<3;i++){
    extr1n(i,3)=t1n[i];
    extr2n(i,3)=t2n[i];
  }
  //instead could just change extrinsics....
  Pn1=newintr*extr1n;
  Pn2=newintr*extr2n;
  
  T1 = Pn1.getCols(first3,3)* (Po1.getCols(first3,3).inverse());
  T2 = Pn2.getCols(first3,3)* (Po2.getCols(first3,3).inverse());
}


//Instantiate some things.
template Image<float> applyHomography(Image<float> ,Matrix &,
				double,double,double,double);
template Image8 applyHomography(Image8,Matrix &,
				double,double,double,double);




