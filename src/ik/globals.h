#ifndef MY_GLOBALS_H
#define MY_GLOBALS_H

#include <nmath/matrix.h>
#include <assert.h>
#include <nmath/vec3.h>
#include <nmath/vec4.h>
#include <nmath/vec2.h>
#include <nimage/image.h>
#include <GL/gl.h>
#include <nmath/mat3x3.h>
#include <nmath/mat4x4.h>

using namespace nacb;
using NIMAGE_NSPACE::Image;

bool readCLB(const char * fname,
	     Mat3x3 & K,
	     Mat4x4 & E);

/** Set the GL_PROJECTION matrix from camera intrinsics.
  Create an openGL projection matrix with the given focal length
  and principle point.  This assumes rendering to an image/viewport
  of size imw,imh and depth range minz,maxz.

  If you write out the equations (assuming viewport of ww,wh):
  x-coord:
  (A11*X/Z+A12*Y/Z+A13*Z/Z+A14+1)*ww/2=a11*(X/Z)+a12*(Y/Z)+a13(Z/Z)
  Assuming A14 is zero, then you get:
  A11=a11*2/ww
  A12=a12*2/ww
  A13=a13*2/ww-1
  (similarily for y)

  The above equation ensures that a point (X,Y,Z) projects to the
  same coordinate using openGL as the typical K*eye(3,4) matrix would.
  
  If viewport is different, say (ww/2)*sc, then it should be scaled properly.
  
  \param fx the (0,0) entry of camera intrinsics
  \param fy the (1,1) entry of the camera intrinsics
  \param px the principle point x-coord(entry (0,2) of camera intrinsics)
  \param py the principle point y-coord(entry (1,2) of camera intrinsics)
  \param imw the image width (the image dimension of the true image, not
         the dimensions of your window)
  \param imh the image height (again the dimension of the camera sensor,
         not the dimensions of your window)
*/
void myProjection(double fx,double fy,
                  double px,double py,
                  int imw,int imh,
                  double minz,double maxz,
		  int origintl=0);

nacb::Imagef unprojectDepth(nacb::Imagef & z, double minz, double maxz, bool origintl = false);

void drawArrow(const nacb::Vec4d & from,const nacb::Vec4d & to,double aspect=0.5,double head=0.9);


GLuint loadProgram(GLuint & vertShader,GLuint & fragShader,
		   const char * vertName="",
		   const char * fragName="");

template <class T>
bool operator==(const nacb::Image<T> & im1,const nacb::Image<T> & im2){
  return (im1.data==im2.data);
}


#endif
