/*
  Load in an image and mesh and compute the disparity.
  
  Brute force, just searches through all the triangles.
*/
#include <algorithm>
#include <vector>

#include "ik/mesh.h"
#include "ik/triangle.h"

#include "autorecon/clbfile.h"
#include "GroundTruthLib.h"
#include "DisparitySequence.h"
#include <nmisc/commandline.h>

using namespace nacb;




int main(int ac, char * av[]){
  nacb::CommandLine cline;
  std::string geomName;
  std::string flowName;
  std::string dispName;
  std::string depthName;
  int geomIndex;

  cline.registerOption("geom", "Geometry base with %%d in it.", &geomName, 'g');
  cline.registerOption("disp", "Output disp file name %%d in it.", &dispName, 'd');
  cline.registerOption("depth", "Output depth file name %%d in it.", &depthName,  0);
  cline.registerOption("flow", "Output flow file name %%d in it.", &flowName, 'f');
  cline.registerOption("gindex", "Geometry index", &geomIndex, 'i');
  cline.parse(ac, av);

  if(optind == ac){
    printf("Need sequence.\n");
    return 0;
  }

  SimpleSequence<unsigned char> seq(av[optind]);

  Mesh mesh;
  mesh.loadobj((boost::format(geomName) % geomIndex).str().c_str());

  std::vector<Triangle> triangles = getTriangles(mesh);
  nacb::Imagef weights;
  nacb::Image<uint32_t> intersects;
  nacb::Imagef depth = computeImageDisparity(seq.image, seq.P, triangles, intersects, weights);

  
  depth.getNormalizedImage().write("/tmp/_gt.png");
  depth.write(dispName.c_str());

  if(depthName.size()){
    nacb::Imagef d = depth.copy();
    for(int i=0; i<depth.w*depth.h; i++){
      d.data[i] = 1.0/d.data[i];
    }
    d.write(depthName.c_str());
  }
    

  if(flowName.size()){
    int index = geomIndex + 1;
    Mesh mesh;

    mesh.loadobj((boost::format(geomName) % index).str().c_str());

    if(mesh.tris.size() != triangles.size()){
      index = geomIndex - 1;
      mesh.loadobj((boost::format(geomName) % index).str().c_str());      
    }
    nacb::Imagef flow = computeImageFlow(seq.image, triangles, getTriangles(mesh),
					 intersects, weights);
        
    if(index < geomIndex)
      flow *= -1;
    
    flow.write(flowName.c_str());
  }

  return 0;
}
