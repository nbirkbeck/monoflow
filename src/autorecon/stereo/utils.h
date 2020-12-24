#ifndef STEREO_UTILS_H
#define STEREO_UTILS_H

#include <GL/gl.h>
#include <vector>
#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include <nmath/matrix.h>
#include <nimage/image.h>
using namespace nacb;

/**
   Note:
   
   Pretty much everything in this file should be merged into some real place.
   Most of them are duplicates with slightly different arguments, or just
   straight up copies from other places.

   One exception is applyHomographyGL...it is new
   
   disparityBounds
**/

//Vec3f backProject(const Matrix & KRinv, const Matrix & Kt, float x,float y,float z);

void rectify(Matrix &, Matrix &,
	     Matrix &, Matrix &,
	     Matrix & T1, Matrix & T2,
	     Matrix & Pn1, Matrix & Pn2);
void findBounds(Image8 & im1,Matrix & H1,
		double & minx,double & maxx);
void findBounds(Image8 & im1,Matrix & H1,
		Image8 & im2,Matrix & H2,
		double & minx,double & maxx,
		double & miny,double & maxy);

void applyHomographyGL(GLuint prog,
		       const Matrix & K,
		       const Matrix & Kinv,
		       const Matrix & dd,
		       const Matrix & H,
		       int w, int h,
		       double minx,
		       double maxx,
		       double miny,
		       double maxy);
void applyHomographyGL(const Matrix & H,
		       double minx,
		       double maxx,
		       double miny,
		       double maxy);

template <class T>
Image<T> applyHomography(Image<T> im,
			 Matrix & H,
			 double minx,
			 double maxx,
			 double miny,
			 double maxy);

nacb::Vec2d getDisparityRange(const nacb::Matrix & Pn1,
			      const nacb::Matrix & Pn2,
			      double minx1, double minx2,
			      double zmin, double zmax);

int readobj_and_cache(const char * fname, 
		      std::vector<Vec3<int> > & tris, 
		      std::vector<Vec3f> & vert, 
		      const nacb::Vec3f & sc = nacb::Vec3f(1,1,1), 
		      const nacb::Vec3f & tr = nacb::Vec3f(0,0,0));

int readobj(const char * fname, 
	    std::vector<Vec3<int> > & tris,
	    std::vector<nacb::Vec3f> & vert,
	    const nacb::Vec3f & sc = nacb::Vec3f(1,1,1), 
	    const nacb::Vec3f & tr = nacb::Vec3f(0,0,0));

int writeobj(const char * fname, std::vector<Vec3<int> > & tris, std::vector<nacb::Vec3f> & vert);



#endif
