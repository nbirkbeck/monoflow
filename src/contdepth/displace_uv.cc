#include "VariationalProblemUV.h"
#include "DisplaceUV.h"
#include "DisparitySequence.h"

#include "utigl/glwindow.h"
#include "ik/mesh.h"
#include "ik/triangle.h"
#include <nimage/image.h>

#include <nmisc/commandline.h>

using namespace nacb;
using namespace std;


nacb::Imagef warpBack(const nacb::Imagef & image, DisplaceUV & warp, const nacb::Imagef & D){
  nacb::Imagef warped(D.w, D.h, image.nchannels);
  for (int y=0; y<D.h; y++) {
    for (int x=0; x<D.w; x++) {
      Vec2f co = warp(x, y, D(x, y));
      
      for (int k=0; k<image.nchannels; k++) {
	warped(x, y, k) = image.bilinear(co.x, co.y, k);
      }
    }
  }
  return warped;
}


int main(int ac, char * av[]){
  std::string geomFile;
  std::string output;

  nacb::CommandLine cline;
  double imageCostWeight = 1.0;
  int downsample = 0;
  int downsampleMax = 0;
  int maxits = 4;
  int medianFilterSize = 0;

  cline.registerOption("geom", "The geometry file", &geomFile, 'g');
  cline.registerOption("img", "Disparity image weight", &imageCostWeight, 'i');
  cline.registerOption("downsample", "Downsample value.", &downsample, 'd');
  cline.registerOption("vdisp_dsmax", "Downsample max.", &downsampleMax, 0);
  cline.registerOption("maxits", "maximum iterations.", &maxits, 0);
  cline.registerOption("output", "Output displacement file.", &output, 0);
  cline.registerOption("medianFilter", "Median filter size.", &medianFilterSize, 0);
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

  printf("%ld %ld %ld\n",
         mesh.vert.size(), mesh.tris.size(), mesh.tvert.size());

  std::vector<nacb::Imagef>  images;
  std::vector<DisplaceUV>     warps;
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

  nacb::Imagef depth(texw, texh, 1);
  depth = 0;

  DisplacedMesh dispMesh;
  BaseGeometry::ptr baseGeometry = createBaseGeometryFromMesh<ImageBaseGeometry>(mesh, texw, texh, dispMesh);

  // Initialize the pairing of images.
  for (int i=0; i<(int)seqs.size(); i++) {
    for (int j=i+1; j<(int)seqs.size(); j++) {
      pairs.push_back(make_pair(i, j));
    }
  }

  int downsample_local = downsample;
  while (downsample_local >= downsampleMax) {    
    BaseGeometry::ptr baseGeometrySubdiv = 
      createBaseGeometryFromMesh<ImageBaseGeometry>(mesh, 
						    texw >> downsample_local, 
						    texh >> downsample_local, dispMesh);
    depth = depth.resize(baseGeometrySubdiv->getWidth(),
			 baseGeometrySubdiv->getHeight());

    // Initialize the warps and the image list.
    images.clear();
    warps.clear();
    weights.clear();

    for (int i=0; i<(int)seqs.size(); i++) {
      images.push_back(seqs[i].image);
      warps.push_back(DisplaceUV(baseGeometrySubdiv, seqs[i].A, seqs[i].E));
      
      weights.push_back(getWeightImage(seqs[i], baseGeometrySubdiv));
    }

    // Solve the problem
    double imageCostWeightDS = imageCostWeight * sqrt(double(1 << downsample_local)) / images.size();

    for (int k=0; k<maxits; k++) {
      if (medianFilterSize)
	depth = depth.medianFilter(medianFilterSize);

      VariationalProblemUV vari(images, warps, weights, pairs, depth, imageCostWeightDS);
      
      nacb::Imagef residuals;
      mg::wcycle_solver solver(3, 3, 10, 10, true);
      
      printf("%lf %lf\n", 
	     VariationalProblemUV::imageDataEnergy(depth, images, warps, weights, pairs),
	     VariationalProblemUV::smoothnessEnergy(depth));
      
      vari.iteration(residuals, false).print();
      solver.solve(vari);
      vari.iteration(residuals, false).print();   
      
      depth = vari.getCombinedResult();
      depth.getNormalizedImage().write("/tmp/res.png");
      
      printf("%lf %lf\n", 
	     VariationalProblemUV::imageDataEnergy(depth, images, warps, weights, pairs),
	     VariationalProblemUV::smoothnessEnergy(depth));
    }
    downsample_local--;
    if (downsample_local > 0)
      seqs = DisparitySequence::downsample(seqs_full_res, downsample_local);
    else if (downsample_local == 0)
      seqs = seqs_full_res;
  }

  // Write the images out.
  for (int k=0; k<(int)seqs.size(); k++) {
    std::string name = (boost::format("/tmp/res-%04d.png") % k).str();
    warpBack(images[k], warps[k], depth).write(name.c_str());
  }
  dispMesh.displace(baseGeometry, depth);
  dispMesh.saveobj("/tmp/displaced.obj", false, true);

  if (output.size())
    depth.write(output.c_str());

  delete window;
  return 0;
}
