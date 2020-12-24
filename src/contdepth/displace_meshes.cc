/**
   A quick program to displace and offset a mesh.

   Neil Birkbeck, 2010
*/

#include "VariationalOffsetProblemUV.h"
#include "DisplaceUV.h"
#include "DisparitySequence.h"

#include "utigl/glwindow.h"
#include "ik/mesh.h"
#include "ik/triangle.h"
#include <nimage/image.h>

#include <nmisc/commandline.h>

using namespace nacb;
using namespace std;

void temporal_box_smooth(std::vector<nacb::Imagef> & images, int maxits = 3) {
  for (int its=0; its<maxits; its++) {
    std::vector<nacb::Imagef> images_out(images.size());

    for (int i=0; i<images.size(); i++) {
      nacb::Imagef image = images[i].copy();
      float sum = 1;

      for (int j=-2; j<=2; j++) {
	if (i + j > 0  && i + j < images.size()) {
	  image += images[i + j];
	  sum++;
	}
      }      
      images_out[i] = image * (1.0/sum);
      images_out[i] = images_out[i].boxFilter(2);
    }
    images = images_out;
  }
}


int main(int ac, char * av[]){
  std::string geomFile;
  std::string dispFile;
  std::string outputFile;

  nacb::CommandLine cline;
  int geomIsBase = 0;
  int base = 0;
  int tsmooth = 0;

  cline.registerOption("geom", "The geometry file", &geomFile, 'g');
  cline.registerOption("disp", "Initialize to this displacement.", &dispFile, 0);
  cline.registerOption("output", "Output name for the objects (should have %d) in it", &outputFile, 0);
  cline.registerOption("geomIsBase", "Geometry should be used as base", &geomIsBase, 0);
  cline.registerOption("tsmooth", "Number of temporal smoothing iterations.", &tsmooth, 0);
  cline.registerOption("base", "Base index used when writing out.", &base, 0);
  cline.parse(ac, av);

  if (!geomFile.size()) {
    std::cout << "Need a geometry file.\n";
    return 2;
  }

  std::cout << "Exporting to output:" << outputFile << "\n";

  Mesh mesh;
  mesh.loadobj(geomFile.c_str());

  printf("%d %d %d\n", mesh.vert.size(), mesh.tris.size(), mesh.tvert.size());

  BaseGeometry::ptr baseGeometry;
  DisplacedMesh dispMesh;
  GLWindow * window = 0;

  std::vector<nacb::Imagef> depths;

  for (int i=optind; i<ac; i++) {
    nacb::Imagef dispFlow(av[i]);
    
    std::cout << "Loading: " << av[i] << "\n";
    std::cout << "Flow: " << dispFlow.w << "," << dispFlow.h << "\n";

    depths.push_back(dispFlow);
  }

  if (tsmooth) temporal_box_smooth(depths, tsmooth);
  
  for (int i=0; i<depths.size(); i++) {
    nacb::Imagef dispFlow = depths[i];

    if (dispFlow.nchannels == 1) {
      nacb::Imagef all(dispFlow.w, dispFlow.h, 4);
      all = 0;
      all.setChannel(0, dispFlow);
      dispFlow = 0;
    }

    if (!baseGeometry.get()) {
      window = new GLWindow(dispFlow.w, dispFlow.h);
      baseGeometry = createBaseGeometryFromMesh(mesh, dispFlow.w, dispFlow.h, dispMesh);
      if (geomIsBase) {
	dispMesh.loadobj(geomFile.c_str());
	dispMesh.loadSourcesFromTextureCoords(dispFlow.w, dispFlow.h);
      }
    }

    dispMesh.displaceAndOffset(baseGeometry, dispFlow);
    dispMesh.saveobj((boost::format(outputFile) % (i + base)).str().c_str(), false, true);

    printf("%ld %ld %ld\n", dispMesh.vert.size(), dispMesh.tris.size(), dispMesh.tvert.size());
  }
  return 0;
}
