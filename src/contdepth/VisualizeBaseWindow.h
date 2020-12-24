#include "AnimationCurve.h"

class VisualizeBaseWindow: public GLWindow {
public:
   VisualizeBaseWindow() : GLWindow(640, 480){    
    font = FFont("Vera.ttf", 12);
    font.scaleToBe(12, height());

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    lighting = true;

    reading = 0;
    readTime = 0;
    
    writingToArchive = 0;
    writeTime = 0;

    screenshot = 0;

    glClearColor(0, 0, 0, 0);
  }

  void doReading() {
    if (reading) {
      if (writeName.size()) {
	nacb::Image8 image(width(), height(), 3);
	glReadPixels(0, 0, width(), height(), image.channelToGL(), image.typeToGL(), image.data);
	image = image.flip();
	image.write((boost::format(writeName) % writeFrameIndex).str().c_str());

	writeFrameIndex++;
      }

      readTime++;

      if (readTime >= curves.max()) {
	reading = false;
	printf("Reading complete: %f\n", curves.max());

	if (writeName.size())
	  exit(1);
      }
      else {
	load(++readTime);
	usleep(30000);
	refresh();
      }
    }
  }

  void doWriting() {
    if (writingToArchive) {
      store(writeTime++);
      usleep(30000);
      refresh();
    }
  }

  void doScreenshot() {
    if (screenshot) {
      static int screenNum = 0;
      std::string screenName = "/tmp/screen-%04d.png";
      nacb::Image8 image(width(), height(), 3);
      glReadPixels(0, 0, width(), height(), 
		   image.channelToGL(), image.typeToGL(), image.data);

      image = image.flip();
      image.write((boost::format(screenName) % screenNum).str().c_str());

      screenNum++;
      screenshot = 0;
    }
  }

  virtual void storeCurveData(double where) {
    curves["lighting"].insert(where, lighting);
    curves["t"].insert(where, t);
    curves["texture"].insert(where, drawTex);

    curves["cpos.x"].insert(where, cpos.x);
    curves["cpos.y"].insert(where, cpos.y);
    curves["cpos.z"].insert(where, cpos.z);

    curves["cquat.x"].insert(where, cquat.v.x);
    curves["cquat.y"].insert(where, cquat.v.y);
    curves["cquat.z"].insert(where, cquat.v.z);
    curves["cquat.a"].insert(where, cquat.a);
  }
  
  virtual void store(double where) {
    // open the archive
    storeCurveData(where);

    std::ofstream ofs("archive.txt");
    boost::archive::text_oarchive oa(ofs);
    oa << (const AnimationCurves &)curves;

    // oa << curves;
  }
  
  virtual bool load(double where){
    if (curves.count("texture"))
      drawTex = (bool) curves["texture"](where);

    if (curves.count("lighting"))
      lighting = (bool) curves["lighting"](where);
    
    if (curves.count("cpos.x"))
      cpos.x = curves["cpos.x"](where);
   
    if (curves.count("cpos.y"))
      cpos.y = curves["cpos.y"](where);

    if (curves.count("cpos.z"))
      cpos.z = curves["cpos.z"](where);

    // Load in the curves for the camera quaternion
    if (curves.count("cquat.x"))
      cquat.v.x = curves["cquat.x"](where);
   
    if (curves.count("cquat.y"))
      cquat.v.y = curves["cquat.y"](where);

    if (curves.count("cquat.z"))
      cquat.v.z = curves["cquat.z"](where);

    if (curves.count("cquat.a"))
      cquat.a = curves["cquat.a"](where);

    if (curves.count("t"))
      t = curves["t"](where);

    cquat.normalize();
    return true;
  }

  void loadArchive(const std::string & name) {
    std::ifstream ifs(name.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> curves;
  }

  void setReading(bool r) {
    reading = r;
    
    if (reading) {
      printf("reading.\n");
      readTime = 0;
      load(readTime);
    }
  }

  bool keyboard(unsigned char c, int x, int y){
    switch (c) {
    case 'R':
      loadArchive("archive.txt");
      setReading(!reading);
      break;

    case 'W':
      writingToArchive = !writingToArchive;
      break;

    case 'w':
      {
	static int wire = 0;
	wire = !wire;
	glPolygonMode(GL_FRONT_AND_BACK, wire ? GL_LINE : GL_FILL);
      }
      break;

    case 'S':
      screenshot = true;
      break;

    case 'l':
      lighting = !lighting;
      break;

    case 't':
      drawTex = !drawTex;
      break;
    }

    refresh();
    return true;
  }

  void setWritingFrames(const std::string & writeName) {
    this->writeName = writeName;
    this->writeFrameIndex = 0;
  }

protected:
  FFont font;

  bool screenshot;
  bool drawTex;
  bool lighting;
  double t;

  AnimationCurves curves;
  double writeTime, readTime;
  bool writingToArchive, reading;

  std::string writeName;
  int writeFrameIndex;
};
