/**
   Quick hack to see what the flow looks like on the object...it ain't good.

   Neil Birkbeck, 2009
 */
#include "utigl/glwindow.h"
#include "utigl/glcommon.h"
#include "utigl/ffont.h"
#include <nmisc/commandline.h>

#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"
#include "autorecon/clbfile.h"

#include <vector>
#include <map>
#include <fstream>

#include <boost/format.hpp>

#include "AnimationCurve.h"
#include "ExplicitDepthOffsetMesh.h"
#include "VisualizeBaseWindow.h"


class Win: public VisualizeBaseWindow {
public:
  Win(const char * imageName,
      const char * calibName,
      const char * depthName,
      const char * flowName, int n = 8){
    printf("Flow: %s\n", flowName);
    meshes = ExplicitDepthOffsetMesh::load(imageName, calibName, depthName, flowName, n);
        
    skip = 4;

    drawMesh = meshes.size();
    showOffs = false;
    drawTex = false;
    backward = false;
    lighting = true;
    drawColors = false;    

    glClearColor(0, 0, 0, 0);
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

    if (drawTex)
      glEnable(GL_TEXTURE_2D);
    else
      glDisable(GL_TEXTURE_2D);

    if (lighting)
      glEnable(GL_LIGHTING);
    else
      glDisable(GL_LIGHTING);

    
    if (animating) {
      ExplicitDepthOffsetMesh::drawBlendedMeshes(meshes, t, skip, showOffs, backward, drawColors);
    }
    else {
      if (drawMesh == (int)meshes.size())
	for(int k=0; k<(int)meshes.size(); k++)
	  meshes[k].draw(skip, showOffs, drawColors);
      else
	meshes[drawMesh].draw(skip, showOffs, drawColors);
    }
    
    doWriting();
    doReading();
    doScreenshot();
  }

  void storeCurveData(double where) {
    VisualizeBaseWindow::storeCurveData(where);

    curves["skip"].insert(where, skip);
    curves["animating"].insert(where, animating);
    curves["showOffs"].insert(where, showOffs);
    curves["drawMesh"].insert(where, drawMesh);
  }
  
  bool load(double where){
    VisualizeBaseWindow::load(where);

    if (curves.count("drawMesh"))
      drawMesh = (int)curves["drawMesh"](where);

    if (curves.count("skip"))
      skip = (int)curves["skip"](where);

    if (curves.count("animating"))
      animating = (bool) curves["animating"](where);

    if (curves.count("showOffs"))
      showOffs = curves["showOffs"](where);
    return true;
  }

  void motion(int x, int y){
    if(bdown == 1){
      t = std::max((double)0, std::min((double)meshes.size()-1, double(y)/height()*(meshes.size() - 1)));
    }
    GLWindow::motion(x, y);
  }
  
  bool keyboard(unsigned char c, int x, int y){
    VisualizeBaseWindow::keyboard(c, x, y);

    switch (c) {
    case 'C':
      {
	nacb::Mat4x4 m = (Matrix::rotx(M_PI)*meshes[0].E).inverse();
	cquat = nacb::Quaternion::fromMatrix(m);
	cpos = nacb::Vec3d(m(0, 3), m(1, 3), m(2, 3));
      }
      break;

    case 'b':
      backward = !backward;
      break;

    case 'c':
      drawColors = !drawColors;
      break;

    case '1':
      skip = 1;
      break;

    case '2':
      skip = 2;
      break;

    case '3':
      skip = 4;
      break;

    case '4':
      skip = 8;
      break;

    case 'S':
      screenshot = true;
      break;

    case '-':
      for(int k=0; k<(int)meshes.size(); k++)
	meshes[k].scaleOffs *= -1;
      break;

    case ' ':
      showOffs = !showOffs;
      break;

    case ']':
      for(int k=0; k<(int)meshes.size(); k++)
	meshes[k].scaleOffs *= 1.2;
      break;
	
    case '[':
      for(int k=0; k<(int)meshes.size(); k++)
	meshes[k].scaleOffs /= 1.2;
      break;

    case '0':
      drawMesh = (drawMesh + 1) % (meshes.size() + 1);
      break;

    case 'a':
      animating = !animating;
      break;
    }

    refresh();
    return true;
  }

protected:
  vector<ExplicitDepthOffsetMesh> meshes;

  bool drawColors;
  int drawMesh;
  bool backward;
  bool animating;
  bool showOffs;
  int skip;
};



int main(int ac, char * av[]){
  if (ac <= 3) {
    printf("Need input image, calibration file, depth and flow\n");
    return 0;
  }

  nacb::CommandLine cline;
  std::string writeName, archiveName;
  int n = 8;
  cline.registerOption("write", "Write the images to this file.", &writeName, 0);
  cline.registerOption("archive", "The archive to read from.", &archiveName, 0);
  cline.registerOption("n", "The number of images to load", &n, 0);
  cline.parse(ac, av);
  

  const char * flowName = av[optind + 3];
  if(std::string("null") == flowName)
    flowName = 0;
  Win * win = new Win(av[optind], av[optind+1], av[optind+2], flowName, n);

  if (cline.gotArgument("write") && cline.gotArgument("archive")) { 
    win->loadArchive(archiveName);
    win->setReading(true);
    win->setWritingFrames(writeName);
  }

  win->loop(1);

  return 0;
}
