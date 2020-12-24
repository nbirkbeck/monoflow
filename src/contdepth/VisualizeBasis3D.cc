/**
   Quick visualization to see the basis flow in 3D.

   Neil Birkbeck, Jan 2011
*/
#include "utigl/glwindow.h"
#include "utigl/glcommon.h"
#include "utigl/ffont.h"
#include <nmisc/commandline.h>
#include "ik/mesh.h"
#include "ik/armature.h"

#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"
#include "autorecon/clbfile.h"

#include <vector>
#include <map>
#include <fstream>

#include <boost/format.hpp>

#include "FlowBasis3D.h"
#include "DisparitySequence.h"
#include "AnimationCurve.h"
#include "VisualizeBaseWindow.h"
#include "BaseMeshAnimation.h"
#include "DisplaceUV.h"
#include "TimeVaryingDisplacedMesh.h"


class Win: public VisualizeBaseWindow {
public:
  Win(char ** av, int ac) {
    seqs = DisparitySequence::loadSequences(av, ac);
    drawTex = false;
    lighting = true;
    glClearColor(0, 0, 0, 0);
    tex = 0;
    activeCamera = 0;
    mesh.drawFlow = 0;
    mesh.drawTangents = 0;

    glGenTextures(1, &cameraTex);
  }

  void drawHud(){    
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0.3);
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-aspect(), aspect(), -1, 1, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glEnable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    double x = -aspect() + 0.05f;
    double y = 0.9f;
    double sp = font.getScaledMaxHeight();
    font.setColor(1, 1, 1, 1);

    char str[1024];
    snprintf(str, sizeof(str), "Time: %f", t);
    font.drawString(str, x, y); 
    y -= sp;

    if (drawTex) {
      font.drawString("Texture", x, y); 
      y -= sp;
    }

    if (lighting) {
      font.drawString("Lighting", x, y);
      y -= sp;
    }

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    
    glPopAttrib();
  }

  void drawScene(){
    drawHud();

    if (seqs.size()) {
      // Drawing the camera image by backprojection ensures that the geometry
      glDisable(GL_LIGHTING);

      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, cameraTex);

      seqs[activeCamera].load(round(t));
      seqs[activeCamera].image.initTexture();
      
      glDepthMask(0);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      //glLoadIdentity();
      //glOrtho(-1, 1, -1, 1, -1, 1);

      nacb::Matrix KRinv = seqs[activeCamera].P.submatrix(0, 0, 3, 3).inverse();
      nacb::Matrix KRt = seqs[activeCamera].P.submatrix(0, 3, 3, 1);

      const float depth = 10.0f;
      float w = seqs[activeCamera].image.w;
      float h = seqs[activeCamera].image.h;

      Vec3f p1 = backProject(KRinv, KRt, 0.f, 0.f, depth);
      Vec3f p2 = backProject(KRinv, KRt, w, 0.f, depth);
      Vec3f p3 = backProject(KRinv, KRt, w, h, depth);
      Vec3f p4 = backProject(KRinv, KRt, 0.f, h, depth);


      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      //glLoadIdentity();

      glColor3f(1, 1, 1);

      glBegin(GL_QUADS);
      glTexCoord2f(0, 0); glVertex3fv(p1.data);
      glTexCoord2f(1, 0); glVertex3fv(p2.data);
      glTexCoord2f(1, 1); glVertex3fv(p3.data);
      glTexCoord2f(0, 1); glVertex3fv(p4.data);
      glEnd();

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();

      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
      glDepthMask(1);
    }

    if (drawTex) {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, tex);
    }
    else
      glDisable(GL_TEXTURE_2D);

    if (lighting)
      glEnable(GL_LIGHTING);
    else
      glDisable(GL_LIGHTING);

    glColor3f(1, 1, 1);
    if (drawMesh)
      mesh.draw();

    doWriting();
    doReading();
    doScreenshot();
    
    mesh.setTime(t);
  }

  void storeCurveData(double where) {
    VisualizeBaseWindow::storeCurveData(where);

    // Store any other curves.
    curves["activeCamera"].insert(where, activeCamera);
    curves["drawMesh"].insert(where, drawMesh);
  }
  
  bool load(double where){
    VisualizeBaseWindow::load(where);

    // Load any other 
    if (curves.count("activeCamera"))
      activeCamera = (int)curves["activeCamera"](where);

    if (curves.count("drawMesh")) {
      drawMesh = (int)curves["drawMesh"](where);
    }
    return true;
  }

  void motion(int x, int y){
    if(bdown == 1){
      t = std::max((double)0, std::min((double)numTimes, double(y)/height()*(numTimes - 1)));
      mesh.setTime(t);
      refresh();
    }
    GLWindow::motion(x, y);
  }

  void setCameraView() {
    if (seqs.size()) {
      nacb::Mat4x4 m = (Matrix::rotx(M_PI)*seqs[activeCamera].E).inverse();
      cquat = nacb::Quaternion::fromMatrix(m);
      cpos = nacb::Vec3d(m(0, 3), m(1, 3), m(2, 3));
    }
    else
      printf("seqs: %d\n", (int)seqs.size());
  }

  bool keyboard(unsigned char c, int x, int y){
    VisualizeBaseWindow::keyboard(c, x, y);

    switch (c) {
    case 'B':
      mesh.drawBaseMesh = !mesh.drawBaseMesh;
      break;

    case 'b':
      mesh.drawBones = !mesh.drawBones;
      break;

    case 'T':
      mesh.drawTangents = !mesh.drawTangents;
      break;

    case 'f':
      mesh.drawFlow = !mesh.drawFlow;
      break;

    case 'm':
      drawMesh = !drawMesh;
      break;

    case '1':
    case '2':
    case '3':
    case '4':
      activeCamera = (c - '1') % seqs.size();
      setCameraView();
      break;

    case 'C':
      {
	setCameraView();
      }
      break;
    }

    refresh();
    return true;
  }

  void setTexture(const std::string& texName) {
    nacb::Image8 texImage;
    texImage.read(texName.c_str());

    if (!tex) {
      glGenTextures(1, &tex);
    }
    glBindTexture(GL_TEXTURE_2D, tex);
    texImage.initTexture();
  }

  int activeCamera;
  TimeVaryingDisplacedMesh mesh;
  int numTimes;
  GLuint tex, cameraTex;
  
protected:
  int drawMesh;
  std::vector<DisparitySequence> seqs;
};



int main(int ac, char * av[]){
  nacb::CommandLine cline;
  int numTimes = 10;
  std::string writeName, archiveName;
  std::string basisFile;
  std::string geomFile, bonesFile, animFile;
  std::string textureFile;

  int texw = 128, texh = 128;

  cline.registerOption("write", "Write the images to this file.", &writeName, 0);
  cline.registerOption("archive", "The archive to read from.", &archiveName, 0);
  cline.registerOption("basisFile", "The basis file.", &basisFile);
  cline.registerOption("geom", "The geometry.", &geomFile);
  cline.registerOption("tex", "The texture.", &textureFile);
  cline.registerOption("bones", "The bones.", &bonesFile);
  cline.registerOption("anim", "The bones.", &animFile);
  cline.registerOption("texw", "The texture width.", &texw);
  cline.registerOption("texh", "The texture height.", &texh);
  cline.registerOption("numTimes", "The number of time frames.", &numTimes);
  cline.parse(ac, av);
  
  Win * win = new Win(av + optind, ac - optind);
  
  win->mesh.setMaximumTextureResolution(texw, texh);
  win->mesh.loadMesh(geomFile.c_str());
  if (!win->mesh.loadBasis(basisFile.c_str())) {
    std::cout << "Using null basis.\n" << std::endl;
    win->mesh.setBasis(FlowBasis3D::ptr(new NullBasis3D()), false);
  }
  win->mesh.loadKinematics(bonesFile, animFile);
  win->mesh.setTime(0);
  if (textureFile.size())
    win->setTexture(textureFile);

  win->numTimes = numTimes;

  if (cline.gotArgument("write") && cline.gotArgument("archive")) { 
    win->loadArchive(archiveName);
    win->setReading(true);
    win->setWritingFrames(writeName);
  }

  win->loop(1);

  return 0;
}
