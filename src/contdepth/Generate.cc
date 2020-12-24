/*
  Quick program to render a synthetic sequence from a given 
  mesh, anim, and a sequence (for calibration).

  Neil Birkbeck, Dec 2, 2009
*/
#include <GL/gl.h>
#include <nmath/mat4x4.h>
#include "ik/mesh.h"
#include "ik/armature.h"
#include "utigl/glwindow.h"
#include "utigl/fbo.h"
#include "utigl/glcommon.h"
#include "autorecon/stereo/simplesequence.h"
#include <nmisc/commandline.h>

class Win : public GLWindow {
public:
  Win() : GLWindow(400, 300) { }
};


int main(int ac, char * av[]){
  nacb::CommandLine cline;
  std::string meshName, bones;
  std::string anim, texName;
  std::string seqName, outputName;
  int nframes = 0;
  GLuint tex;

  cline.registerOption("mesh", "mesh", &meshName, 0);
  cline.registerOption("bones", "Bones", &bones, 0);
  cline.registerOption("anim", "Animation file", &anim, 0);
  cline.registerOption("nframes", "Number of frames", &nframes, 0);
  cline.registerOption("tex", "Texture name", &texName, 0);
  cline.registerOption("seq", "Sequence name", &seqName, 0);
  cline.registerOption("oimage", "Output imgae name", &outputName, 0);
  cline.parse(ac, av);
  
  Mesh mesh;
  Armature arm;
  PoseKeys keys;

  mesh.loadobj(meshName.c_str());
  arm.read(bones.c_str());
  keys.load(anim.c_str());

  SimpleSequence<unsigned char> seq(seqName);

  Win* win = new Win();
  (void)win;

  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  (nacb::Image8(texName.c_str())).initTexture();

  FrameBufferObject fbo(1024, 1024);
  
  fbo.bind(1);

  for(int i=0; i<nframes; i++){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    seq.load(i);

    keys.setCurrent(i);
    arm.set(keys);
    arm.animate(mesh);

    nacb::Matrix K = seq.A;
    nacb::Matrix E = seq.E;
    E.printMatlab("E");
    
    std::cout << mesh.vert.size() << "\n";

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    cameraPerspectiveMatrix(K(0,0), K(1,1), K(0,2), K(1,2), 
			    seq.image.w, seq.image.h, 0.1, 100, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMultMatrixd(E.transpose().data);

    glBindTexture(GL_TEXTURE_2D, tex);
    glEnable(GL_TEXTURE_2D);

    mesh.draw();

    //fbo.bind(0);

    nacb::Image8 im(seq.image.w, seq.image.h, 4);
    glReadPixels(0, 0, im.w, im.h, im.channelToGL(), im.typeToGL(), im.data);

    if(outputName.size())
      im.write((boost::format(outputName) % i).str().c_str());

    fbo.bind(1);
  }
 
  return 0;
}
