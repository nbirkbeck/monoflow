#include "OffsetMesh.h"
#include "utigl/glcommon.h"

#ifndef glGenBuffers 
#define glGenBuffers glGenBuffersARB
#endif

OffsetMesh::OffsetMesh(const nacb::Image8 & image,
		       const nacb::Imagef & depth,
		       const nacb::Imagef & offset,
		       const nacb::Matrix & P) : DepthMesh(image, depth, P) {
  glGenTextures(1, &offsetTex);
  setOffsetImage(offset);

  // Delete any of the programs created by depth-mesh.

  glDetachShader(prog, fragShader);
  glDetachShader(prog, vertShader);
  glDeleteProgram(prog);
  glDeleteShader(fragShader);
  glDeleteShader(vertShader);

  prog = load_program(vertShader, fragShader, "offset_mesh.vsh", "depth_mesh.fsh");

  offsetScale = 1;
  // setShading(1);
}


void OffsetMesh::setOffsetImage(const nacb::Imagef & offset) {
  glBindTexture(GL_TEXTURE_2D, offsetTex);

  glTexImage2D(GL_TEXTURE_2D,0, GL_RGBA32F_ARB, offset.w, offset.h,
	       0,GL_RGB,GL_FLOAT, offset.data);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);

  glBindTexture(GL_TEXTURE_2D, 0);
}


void OffsetMesh::draw(){
  glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT);

  // Setup the rendering.
  glActiveTexture(GL_TEXTURE0);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, tex);
    
  glActiveTexture(GL_TEXTURE1);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, depthTex);

  glActiveTexture(GL_TEXTURE2);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, offsetTex);

  // Setup the shader parameters.
  glUseProgram(prog);
  
  DepthMesh::setProgramParameters();
  GLint loc = glGetUniformLocation(prog, "offsetScale");
  glUniform1f(loc, offsetScale);

  loc = glGetUniformLocation(prog, "Einv");
  nacb::Matrix Einv = E.inverse();
  float Einvdata[16];
  for(int i=0; i<16; i++)
    Einvdata[i] = Einv.data[i];
  glUniformMatrix4fv(loc, 1, true, Einvdata);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  // DepthMesh::multModelviewMatrix();
  
  if(mesh)
    mesh->draw();

  glPopMatrix();

  glPopAttrib();
}


#ifdef TEST_OFFSET_MESH
/**
   Mostly a clone of test_depth_mesh (it looks like the offsets are 
   correct).
*/

#include "autorecon/clbfile.h"
#include "utigl/glwindow.h"

using namespace nacb;

class Win: public GLWindow {
public:
  Win(const std::string & imageName,
      const std::string & calibName,
      const std::string & depthName,
      const std::string & flowName){
    Matrix K(3, 3), P(3, 4), d(5, 1);
    E = Matrix::eye(4, 4);
    ClbFile::read(calibName.c_str(), K.data, E.data, P.data, d.data);    
    P = K*Matrix::eye(3, 4)*E;
    E.printMatlab("E");

    offsetMesh = new OffsetMesh(Image8(imageName.c_str()),
				Imagef(depthName.c_str()),
				Imagef(flowName.c_str()), P);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER,0.01);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  virtual bool keyboard(unsigned char c, int x, int y){
    switch(c){
    case 27:
      exit(0);
      break;

    case 'c':
      {
	nacb::Mat4x4 m = (Matrix::rotx(M_PI)*E).inverse();
	cquat = nacb::Quaternion::fromMatrix(m);
	cpos = nacb::Vec3d(m(0, 3), m(1, 3), m(2, 3));
      }
      break;

    case '+':
    case '=':
      offsetMesh->setSpacing(offsetMesh->getSpacing()+1);
      break;

    case '_':
    case '-':
      offsetMesh->setSpacing(offsetMesh->getSpacing()-1);
      printf("spacing: %d\n", offsetMesh->getSpacing());
      break;

    case 'a':
      offsetMesh->setShading(1.0 - offsetMesh->getShading());

    case '1':
      offsetMesh->setOffsetScale(-1);
      break;

    case '2':
      offsetMesh->setOffsetScale(0);
      break;

    case '3':
      offsetMesh->setOffsetScale(1);
      break;

    case 'w':
      {
	static int wire = 0;
	wire = !wire;
	if(wire)
	  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
	  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      }
      break;      
    }

    refresh();
    return false;
  }

  ~Win(){
    delete offsetMesh;
  }

  void drawScene(){
    offsetMesh->draw();
  }

  Matrix E;
  OffsetMesh * offsetMesh;
};

int main(int ac, char * av[]){
  if(ac < 5){
    fprintf(stderr, "Need image, calib, depth, and flow.\n");
    return 0;
  }
  Win * win = new Win(av[1], av[2], av[3], av[4]);

  win->loop(1);

  return 0;
}

#endif
