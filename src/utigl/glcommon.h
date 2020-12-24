#ifndef NGL_COMMON_H
#define NGL_COMMON_H

#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES
#endif
#include <GL/gl.h>
#include <sys/types.h>
#include <string>
#include <vector>

std::string resolvePath(const std::string & fname,
		   const std::vector<std::string> & paths);
time_t getModifiedTime(const std::string & fname);

void printShaderInfoLog(GLuint shader,const char * name);
void printProgramInfoLog(GLuint handle,const char * name);

time_t load_shader(const std::string & name,GLuint shader);
void checkLinkStatus(GLuint prog,const char * progname,
		     const char * filename,const char * function,int lineno);
time_t load_arb_program(const std::string & name,GLenum target);

bool link_program(GLuint prog, GLuint s1, GLuint s2, GLuint s3=0, GLuint s4=0);
GLuint create_linked_program(GLuint s1, GLuint s2, GLuint s3=0, GLuint s4=0);

GLuint load_fsh(const std::string & name);
GLuint load_vsh(const std::string & name);

GLuint load_program(GLuint & vertShader,
		    GLuint & fragShader,
		    const std::string & vertName="",
		    const std::string & fragName="");

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
void cameraPerspectiveMatrix(double fx,double fy,
			     double px,double py,
			     int imw,int imh,
			     double minz,double maxz,
			     int origintl=0);



GLenum printGLError_helper(const char * filename,int line);
#define printGLError() printGLError_helper(__FILE__,__LINE__)


#endif
