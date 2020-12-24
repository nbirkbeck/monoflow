#include <GL/glew.h>
#include <GL/gl.h>
#include <nimage/image.h>
#include <nmisc/commandline.h>
#include "autorecon/stereo/simplesequence.h"
#include "ExplicitDepthOffsetMesh.h"
#include <math.h>

#include "utigl/glwindow.h"
#include "utigl/fbo.h"

#include "utils.h"

class Win : public GLWindow {
public:
  Win() : GLWindow(400, 300) {
    t = 0;
    drawImage = false;
    flip = false;
    backward = false;
    finish = false;
    
    fbo = 0;

    glEnable(GL_TEXTURE_2D);    
  }

  void applyModelview() {
    Matrix E = sequence->E;
    if(flip) {
      glRotatef(180, 1, 0, 0);
    }
    glMultMatrixd(E.transpose().data);
    
  }

  void applyProjection() {
    Matrix K = sequence->A;
    int w = sequence->image.w;
    int h = sequence->image.h;
    
    double px = flip ? w - K(0, 2): K(0, 2);
    cameraPerspectiveMatrix(K(0,0), K(1,1), px, K(1,2), w, h, 0.1, 100, !flip);
  }

  void setFrame(int val) {
    t = val;
    if (t >= (int)render.size())
      t = 0;

    /*bool loaded =*/
    sequence->load(render[t].second);    
    // printf("Loaded: %d. %d\n", loaded, render[t].first);
  }

  void drawScene() {
    if(!meshes) 
      return;

    int skip = 1;
    bool showOffs = true;
    bool drawColors = false;

    if (drawImage) {
      static GLuint tex = 0;
      
      if (tex == 0){
	glGenTextures(1, &tex);
      }
      nacb::Image8 & image = sequence->image;
      glBindTexture(GL_TEXTURE_2D, tex);
      image.initTexture();

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      glOrtho(0, 1, 0, 1, -1, 1);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      
      glBindTexture(GL_TEXTURE_2D, tex);

      double top = flip;

      glColor3f(1, 1, 1);
      glBegin(GL_QUADS);
      glTexCoord2f(0, top);
      glVertex2f(0, 0);

      glTexCoord2f(1, top);
      glVertex2f(1, 0);
      
      glTexCoord2f(1, 1 - top);
      glVertex2f(1, 1);
      
      glTexCoord2f(0, 1 - top);
      glVertex2f(0, 1);
      glEnd();      

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();

      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
    }
    else 
      ExplicitDepthOffsetMesh::drawBlendedMeshes(*meshes, render[t].first, skip, 
						 showOffs, backward, drawColors);

    nacb::Image8 results(width(), height(), 4);
    results = 0;

    if(fbo){
      nacb::Image8 resBig(fbo->getWidth(), fbo->getHeight(), 4);
      fbo->bind(0);
      
      glBindTexture(GL_TEXTURE_2D, fbo->getColorTexture());
      glGetTexImage(GL_TEXTURE_2D, 0, resBig.channelToGL(), resBig.typeToGL(), resBig.data);
      glBindTexture(GL_TEXTURE_2D, 0);
      
      results = resBig.subimage(0, 0, width(), height());
      
      fbo->bind(1);
    }
    else 
      glReadPixels(0, 0, width(), height(), results.channelToGL(), results.typeToGL(), results.data);

    if (writeFrameName.size()) {
      results.save((boost::format(writeFrameName) % t).str().c_str());
    }
	

    image_stats_t s;
    std::cout << (s = image_compare(results.getChannel(3), results, sequence->image)) << "\n";;

    if (t >= (int)stats.size())
      stats.push_back(s);

    stats[t] = s;

    if (stats.size() == render.size() && finish) {
      double mean = 0;
      printf("stats=[");
      for(int k=0; k<(int)stats.size(); k++) {
	mean += stats[k].mean;
	printf("%lf,", stats[k].mean);
      }
      printf("];\n\n");
      printf("%f\n", mean/stats.size());
      exit(1);
    }
      

    setFrame(t + 1);

    // usleep(100000);
    refresh();
  }

  bool keyboard(unsigned char c, int x, int y){ 
    switch (c) {
    case 'F':      
      finish = true;
      break;

    case 'f':
      flip = !flip;
      break;

    case 'b':
      backward = !backward;
      break;

    case 'I':
      drawImage = !drawImage;
      break;
    }
    return true;
  }

  FrameBufferObject * fbo;

  int t;

  bool finish;
  bool backward;

  bool flip;
  bool drawImage;

  std::vector<image_stats_t> stats;

  SimpleSequence<unsigned char> * sequence;
  std::vector<ExplicitDepthOffsetMesh> * meshes;
  std::vector<std::pair<double, int> >   render;

  std::string writeFrameName;
};



/*
  Given a directory of *.clb, find those that match those
  "calib.clb" files in the directory style structure.  
*/
int main(int ac, char * av[]){
  nacb::CommandLine cline;
  std::string seqName = "";
  std::string imageName = "";
  std::string flowName = "";
  std::string depthName = "";
  std::string calibName = "";
  std::string writeFrameName = "";
  int visualize = false;
  int backward = false;
  int n = 6;
  int drawImage = 0;
  
  cline.registerOption("backward", "Use forward/backward", &backward, 0);
  cline.registerOption("visualize", "", &visualize, 0);
  cline.registerOption("image", "Image base.", &imageName, 0);
  cline.registerOption("calib", "Calib base name", &calibName, 0);
  cline.registerOption("depth", "Depth base name.", &depthName, 0);
  cline.registerOption("flow", "Flow base name.", &flowName, 0);
  cline.registerOption("drawImage", "Draw the image", &drawImage, 0);
  cline.registerOption("write", "Write out the image frames.", &writeFrameName, 0);

  cline.registerOption("n", "The number of images in the flow/offset sequence to load.", &n, 0);
  cline.registerOption("seq", "The sequnece to load for comparison", &seqName, 0);
  cline.parse(ac, av);

  if (!seqName.size() || !imageName.size() ) {
    cline.printHelp();
    return 0;
  }

  printf("Using sequence: %s\n", seqName.c_str());

  std::cout << "Image: " << imageName << "\n";
  std::cout << "Calib: " << calibName << "\n";
  std::cout << "Depth: " << depthName << "\n";
  std::cout << "Flow : " << flowName << "\n";

  try { 
    Win * window = new Win();
    SimpleSequence<unsigned char> sequence(seqName);

    std::cout << "Loading meshes: " << flowName.c_str() << "\n";
    std::vector<ExplicitDepthOffsetMesh> meshes = 
      ExplicitDepthOffsetMesh::load(imageName.c_str(),
				    calibName.c_str(),
				    depthName.c_str(),
				    flowName.size() ? flowName.c_str() : 0, n);
    std::vector<int> matches;    
    std::vector<double>    t;
    int lower = 0;

    std::cerr << "Got " << n << " meshess.\n";

    for (int i=0; i<(int)meshes.size(); i++) {      
      while (sequence.load(lower)) {
	Matrix diff = sequence.E - meshes[i].E;

	sequence.E.printMatlab("Ese");
	meshes[i].E.printMatlab("Eme");

	if (diff.dot(diff) < 1e-10) {
	  printf("%d maps to %d\n", i, lower);
	  matches.push_back(lower);
	  break;
	}
	else 
	  printf("trying %d to %d..\n", i, lower);
	lower ++;
      }
      if(!sequence.load(lower)) {
	std::cout << "Cannot find matching image.\n";
	break;
      }
    }

    // This is a list of the matching pairs (time to integer index).
    std::vector<std::pair<double, int> > render;

    for (int i=0; i<((int)matches.size() - 1); i++) {
      double t0 = double(i);
      double t1 = double(i + 1);
      
      for (int j=matches[i]; j<matches[i+1]; j++){
	double textra = double(j - matches[i])/(matches[i+1] - matches[i]);
	double t = t0 + textra*(t1 - t0);

	render.push_back(std::pair<double, int>(t, j));

	printf("%f\n", t);
      }
    }

    FrameBufferObject fbo(1024, 1024);
    
    window->sequence = &sequence;
    window->meshes = &meshes;
    window->render = render;
    window->setFrame(0);
    window->backward = backward;

    if(!visualize){      
      fbo.bind(1);
      window->finish = true;
      window->fbo = &fbo;
      window->writeFrameName = writeFrameName;
      window->drawImage = drawImage;
    }

    window->loop(1);
  }
  catch (std::string & e) {
    std::cout << "string exception  " <<  e << "\n";
  }
  catch (exception & e) {
    std::cout << typeid(e).name() << ":" << e.what();
  }
  return 0;
}
