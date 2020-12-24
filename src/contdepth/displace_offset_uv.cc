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


nacb::Imagef warpBack(const std::vector<nacb::Imagef> & images, 
		      DisplaceOffsetUV & warp, const nacb::Imagef & D, int t = 0){
  nacb::Imagef warped(D.w, D.h, images[t].nchannels);
  for (int y=0; y<D.h; y++) {
    for (int x=0; x<D.w; x++) {
      Vec2f co = warp(x, y, D(x, y, 0), Vec3f(D(x, y, 1), D(x, y, 2), D(x, y, 3))*t);
      
      for (int k=0; k<images[t].nchannels; k++) {
	warped(x, y, k) = images[t].bilinear(co.x, co.y, k);
      }
    }
  }
  return warped;
}


int main(int ac, char * av[]){
  std::string geomFile;
  std::string dispFile;
  std::string outputFile;

  nacb::CommandLine cline;
  double imageCostFlowWeight = 1.0;
  double imageCostStereoWeight = 1.0;

  int downsample = 0;
  int downsampleMax = 0;
  int maxits = 4;

  int index1 = 0;
  int index2 = 75;

  cline.registerOption("geom", "The geometry file", &geomFile, 'g');
  cline.registerOption("img", "Disparity image weight (stereo term, will be normalized by pairs)", &imageCostStereoWeight, 'i');
  cline.registerOption("flow", "Flow image weight (will be normalized by number of seqs)", &imageCostFlowWeight, 'f');

  cline.registerOption("downsample", "Downsample value.", &downsample, 'd');
  cline.registerOption("vdisp_dsmax", "Downsample max.", &downsampleMax, 0);
  cline.registerOption("index1", "Image index 1.", &index1, 0);
  cline.registerOption("index2", "Image index 2.", &index2, 0);
  cline.registerOption("maxits", "Maximum iterations.", &maxits, 0);
  cline.registerOption("disp", "Initialize to this displacement.", &dispFile, 0);
  cline.registerOption("output", "Output displacement name.", &outputFile, 0);
  cline.parse(ac, av);

  if (optind + 1 >= ac) {
    std::cout << "Need at least two image sequences.\n";
    return 1;
  }

  if (!geomFile.size()) {
    std::cout << "Need a geometry file.\n";
    return 2;
  }

  Mesh mesh;
  mesh.loadobj(geomFile.c_str());

  printf("%ld %ld %ld\n", mesh.vert.size(), mesh.tris.size(), mesh.tvert.size());

  std::vector<std::vector<nacb::Imagef> >  images;
  std::vector<DisplaceOffsetUV>     warps;
  std::vector<nacb::Imagef>    weights;
  std::vector<pair<int, int> >   pairs;
  
  std::vector<DisparitySequence> seqs_full_res = 
    DisparitySequence::loadSequences(av + optind, ac - optind);
  std::vector<DisparitySequence> seqs;

  if(downsample)
    seqs = DisparitySequence::downsample(seqs_full_res, downsample);
  else
    seqs = seqs_full_res;

  int texw = 512, texh = 512;
  GLWindow * window = new GLWindow(texw, texh);

  nacb::Imagef depthOffset(texw, texh, 4);
  depthOffset = 0;

  if (dispFile.size()) {
    depthOffset = nacb::Imagef(dispFile.c_str());
    if (depthOffset.nchannels <= 1) {
      nacb::Imagef everything(depthOffset.w, depthOffset.h, 4);
      everything = 0;
      everything.setChannel(0, depthOffset);
      depthOffset = everything;
    }
    depthOffset = depthOffset.resize(depthOffset.w >> downsample,
				     depthOffset.h >> downsample);				     
  }

  typedef ImageBaseGeometry BaseGeometryType;
  DisplacedMesh dispMesh;
  BaseGeometry::ptr baseGeometry = createBaseGeometryFromMesh<BaseGeometryType>(mesh, texw, texh, dispMesh);

  // Initialize the pairing of images.
  for (int i=0; i<(int)seqs.size(); i++) {
    for (int j=i+1; j<(int)seqs.size(); j++) {
      printf("Adding pair: %d %d\n", i, j);
      pairs.push_back(make_pair(i, j));
    }
  }

  


  int downsample_local = downsample;
  while (downsample_local >= downsampleMax) {    
    BaseGeometry::ptr baseGeometrySubdiv = 
      createBaseGeometryFromMesh<BaseGeometryType>(mesh, 
				 texw >> downsample_local, 
				 texh >> downsample_local, dispMesh);
    depthOffset = depthOffset.resize(baseGeometrySubdiv->getWidth(),
				     baseGeometrySubdiv->getHeight());

    // Initialize the warps and the image list.
    images.clear();
    warps.clear();
    weights.clear();

    for (int i=0; i<(int)seqs.size(); i++) {
      int w = seqs[i].image.w, h = seqs[i].image.h;
      vector<Imagef> imagesi;
      seqs[i].load(index1);
      seqs[i].image = seqs[i].image.resize(w, h);
      imagesi.push_back(seqs[i].image);
	    
      seqs[i].load(index2);
      seqs[i].image = seqs[i].image.resize(w, h);
      imagesi.push_back(seqs[i].image);

      images.push_back(imagesi);
      warps.push_back(DisplaceOffsetUV(baseGeometrySubdiv, seqs[i].A, seqs[i].E));
      weights.push_back(getWeightImage(seqs[i], baseGeometrySubdiv));

      seqs[i].A.printMatlab("A");

      printf("writing %d\n", i);
      weights.back().write((boost::format("/tmp/weights-%04d.png") % i).str().c_str());
    }

    printf("Solving problem.\n");
    // Solve the problem
    Vec2f imageCostWeightPairDS(imageCostStereoWeight * sqrt(double(1 << downsample_local)) / pairs.size(),
				imageCostFlowWeight * sqrt(double(1 << downsample_local)) / images.size());

    for (int k=0; k<maxits; k++) {
      VariationalOffsetProblemUV vari(images, warps, weights, 
				      pairs, depthOffset, imageCostWeightPairDS);
      printf("Created.\n");
      nacb::Imagef residuals;
      mg::wcycle_solver solver(2, 2, 4, 4, true);
      
      
      printf("%lf %lf \n", 
	     VariationalOffsetProblemUV::imageDataEnergy(depthOffset, images, warps, weights, pairs),
	     VariationalOffsetProblemUV::smoothnessEnergy(depthOffset));
	     
      
      vari.iteration(residuals, false).print();
      solver.solve(vari);
      vari.iteration(residuals, false).print();   
      
      depthOffset = vari.getCombinedResult();

      depthOffset.getChannel(0).getNormalizedImage().write("/tmp/res.png");
      depthOffset.getChannel(1).getNormalizedImage().write("/tmp/res-u.png");
      depthOffset.getChannel(2).getNormalizedImage().write("/tmp/res-v.png");
      depthOffset.getChannel(3).getNormalizedImage().write("/tmp/res-w.png");
      
      // Write the images out.
      for (int k=0; k<(int)seqs.size(); k++) {
	for (int t=0; t<=1; t++) {
	  std::string name = (boost::format("/tmp/res-%04d-%d.png") % k % t).str();
	  warpBack(images[k], warps[k], depthOffset, t).write(name.c_str());
	}
      }
      
      printf("%lf %lf\n", 
	     VariationalOffsetProblemUV::imageDataEnergy(depthOffset, images, warps, weights, pairs),
	     VariationalOffsetProblemUV::smoothnessEnergy(depthOffset));
    }
    downsample_local--;
    if (downsample_local > 0)
      seqs = DisparitySequence::downsample(seqs_full_res, downsample_local);
    else if (downsample_local == 0)
      seqs = seqs_full_res;
  }

  BaseGeometry::ptr baseGeometrySubdiv = 
    createBaseGeometryFromMesh<BaseGeometryType>(mesh, 
						 texw >> downsample_local, 
						 texh >> downsample_local, dispMesh);
  depthOffset = depthOffset.resize(baseGeometrySubdiv->getWidth(),
				   baseGeometrySubdiv->getHeight());
  
  dispMesh.displace(baseGeometrySubdiv, depthOffset);
  dispMesh.saveobj("/tmp/displaced.obj", false, true);

  dispMesh.displaceAndOffset(baseGeometrySubdiv, depthOffset);
  dispMesh.saveobj("/tmp/displaced_offset.obj", false, true);

  if (outputFile.size()) {
    depthOffset.resize(texw, texh).write(outputFile.c_str());
    depthOffset.resize(texw, texh).getChannel(0).getNormalizedImage().write("/tmp/d-final.png");
    depthOffset.resize(texw, texh).getChannel(2).getNormalizedImage().write("/tmp/dv-final.png");
  }

  delete window;
  return 0;
}
