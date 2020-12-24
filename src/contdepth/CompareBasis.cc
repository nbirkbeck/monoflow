/**
   A function to compare the results of the basis flow estimation to the ground truth.

   We want, errors on disparity (mean/median), errors on the flows (mean/median, per frame,
   and averaged over all frames).

   Neil Birkbeck, Jan 2011
*/
#include "utils.h"
#include "TimeVaryingDisplacedMesh.h"
#include "DisplaceOffsetBasis.h"
#include "utigl/glwindow.h"
#include "utigl/fbo.h"
#include <nmisc/commandline.h>
#include <nimage/image.h>
#include <vector>


nacb::Imagef convertPositionToDisplacement(const BaseGeometry::ptr& baseGeometry,
					   const nacb::Imagef& displacementOffset) {
  nacb::Imagef image(displacementOffset.w, displacementOffset.h, 1);

  for (int y=0; y<displacementOffset.h; y++) {
    for (int x=0; x<displacementOffset.w; x++) {
      if (baseGeometry->getMask()(x, y)) {
	nacb::Vec3f normal;
	nacb::Vec3f co = baseGeometry->displace(x, y, 0, &normal);
	
	nacb::Vec3f displaced(displacementOffset(x, y, 0),
			      displacementOffset(x, y, 1),
			      displacementOffset(x, y, 2));
	double depth = normal.dot(displaced - co);
	image(x, y, 0) = depth;
      }
      else 
	image(x, y, 0) = 0;
    }
  }
  return image;
}

/** \brief Take a base geometry and a displacement/offset, and create an image of
           3D coordinates.
 */
nacb::Imagef convertDisplacementToPosition(const BaseGeometry::ptr& baseGeometry,
					   const nacb::Imagef& displacementOffset,
					   nacb::Imagef& worldDisplacement) 
{
  nacb::Imagef image(displacementOffset.w, displacementOffset.h, 3);
  worldDisplacement = nacb::Imagef(image.w, image.h, 3);

  for (int y=0; y<displacementOffset.h; y++) {
    for (int x=0; x<displacementOffset.w; x++) {
      nacb::Vec3f co = baseGeometry->displaceAndOffset(x, y, displacementOffset(x, y, 0),
						       displacementOffset(x, y, 1),
						       displacementOffset(x, y, 2),
						       displacementOffset(x, y, 3));
      image(x, y, 0) = co.x;
      image(x, y, 1) = co.y;
      image(x, y, 2) = co.z;

      nacb::Vec3f rest = baseGeometry->displace(x, y, 0);
      worldDisplacement(x, y, 0) = co.x - rest.x;
      worldDisplacement(x, y, 1) = co.y - rest.y;
      worldDisplacement(x, y, 2) = co.z - rest.z;
    }
  }
  return image;
}


static void imageComparison(std::vector<DisparitySequence> & seqs,
			    TimeVaryingDisplacedMesh & recon, 
			    TimeVaryingDisplacedMesh & truth,
			    int numTimes) {
  std::vector<int> times;
  for (int i=0; i<numTimes; i++)
    times.push_back(i);

  std::vector<std::vector<nacb::Imagef> >  images;
  std::vector<DisplaceOffsetUV>             warps;
  std::vector<std::vector<nacb::Imagef> > weights;

  // Grab all images and initialize the warps.
  DisplaceOffsetBasis::initImagesWarpsAndWeights(seqs, recon.basis, recon.baseMeshAnimation, 
						 recon.dispMesh, recon.baseGeometry,
						 times, recon.basisData, images, warps, weights, true);

  // Warp back all images.
  std::vector<std::vector<nacb::Imagef> > warped;
  for (int i=0; i<(int)seqs.size(); i++) {
    std::vector<nacb::Imagef> warpedThisSeq;

    for (int j=0; j<(int)times.size(); j++) {
      warpedThisSeq.push_back(DisplaceOffsetBasis::warpBack(images[i], recon.basis, recon.baseMeshAnimation, warps[i], recon.basisData, times[j]));
      (warpedThisSeq.back()).write((boost::format("/tmp/comp-%04d-%d.png") % i % j).str().c_str());
    }
    warped.push_back(warpedThisSeq);
  }


  image_stats_t totalFlow;
  image_stats_t totalStereo;

  // Compute the flowed image stats.
  for (int i=0; i<(int)seqs.size(); i++) {
    image_stats_t stats;

    for (int j=0; j<int(warped[i].size()-1); j++) {
      stats += image_compare(weights[i][j] * weights[i][j+1], warped[i][j], warped[i][j+1]);
    }
    stats /= (warped[i].size()-1);
    std::cout << "seq" << i << "(flow):" << stats << std::endl;
    totalFlow += stats;
  }


  // Compute the stereo stats;
  for (int t=0; t<(int)times.size(); t++) {
    image_stats_t stats;
    int count = 0;

    for (int j1=0; j1<(int)seqs.size(); j1++) {
      for (int j2=j1 + 1; j2<(int)seqs.size(); j2++) {
	stats += image_compare(weights[j1][t] * weights[j2][t], warped[j1][t], warped[j2][t]);
	count++;
      }
    }
    stats /= count;

    std::cout << "stereo" << t << ":" << stats << std::endl;
    totalStereo += stats;
  }

  
  // Average the stats.
  totalStereo /= times.size();
  totalFlow /= seqs.size();

  // Output
  std::cout << "total-flow:" << totalFlow << std::endl;
  std::cout << "total-stereo:" << totalStereo << std::endl;
}



int main(int ac, char * av[]){
  nacb::CommandLine cline;
  int numTimes = 10;
  std::string writeName, archiveName;
  std::string basisFile;
  std::string geomFile, bonesFile, animFile;
  std::string truthBasisFile;
  int zeroDisplacements = 0;
  int texw = 512, texh = 512;

  cline.registerOption("basisFile", "The basis file.", &basisFile);
  cline.registerOption("truthBasis", "The basis file for the ground truth.", &truthBasisFile);
  cline.registerOption("geom", "The geometry.", &geomFile);
  cline.registerOption("bones", "The bones.", &bonesFile);
  cline.registerOption("anim", "The bones.", &animFile);
  cline.registerOption("numTimes", "The number of time frames.", &numTimes);
  cline.registerOption("texw", "Texture width", &texw);
  cline.registerOption("texh", "Texture height", &texh);
  cline.registerOption("zero", "Zero the displacements.", &zeroDisplacements, 0);
  cline.parse(ac, av);
  
  GLWindow window(16, 16);
  FrameBufferObject fbo(1024, 1024);
  fbo.bind(1);

  TimeVaryingDisplacedMesh recon, truth;

  recon.loadMesh(geomFile.c_str());
  recon.loadBasis(basisFile.c_str());
  recon.resizeBasisData(texw, texh);
  recon.loadKinematics(bonesFile, animFile);

  truth.loadMesh(geomFile.c_str());
  truth.loadBasis(truthBasisFile.c_str());  
  truth.resizeBasisData(texw, texh);
  truth.loadKinematics(bonesFile, animFile);

  std::cout << recon.baseMeshAnimation.getBaseGeometry().get() << " " << 
    recon.baseGeometry.get() << std::endl;

  //recon.baseMeshAnimation.setBaseGeometry(recon.baseGeometry);
  //truth.baseMeshAnimation.setBaseGeometry(truth.baseGeometry);
  


  nacb::Imagef mask = recon.getBaseGeometry()->getMask();
  nacb::Image8 mask8;
  mask8 = mask;

  if (recon.basisData.w != truth.basisData.w) {
    std::cout << "Flow sizes do not agree, set texw/texh." << std::endl;
    exit(1);
  }

  std::cout << "\n\n\n\n\n";

  std::cout << "First channel comparison.\n";
  // Check the depth error.
  truth.setTime(0);
  recon.setTime(0);

  // First, print out some depth stats (only relevent if first channel is depth).
  nacb::Imagef depthTruth = truth.basis->getDisplacementMap(truth.basisData, 0);
  nacb::Imagef depthRecon = recon.basis->getDisplacementMap(recon.basisData, 0);

  double n, x;
  depthTruth.getChannel(0).getRange(n, x);
  std::cout << "Truth depth range:" << n << " " << x << std::endl;

  depthRecon.getChannel(0).getRange(n, x);
  std::cout << "Recon depth range:" << n << " " << x << std::endl;

  depthTruth.getChannel(0).getNormalizedImage().write("/tmp/t.png");
  depthRecon.getChannel(0).getNormalizedImage().write("/tmp/r.png");

  image_stats_t truth_zero = disp_compare(mask8, truth.basisData.getChannel(0), recon.basisData.getChannel(0)*0);
  std::cout << "truth (zero):" << truth_zero << std::endl;

  image_stats_t recon_zero = disp_compare(mask8, recon.basisData.getChannel(0), truth.basisData.getChannel(0)*0);
  std::cout << "recon (zero):" << recon_zero << std::endl;

  image_stats_t stats = disp_compare(mask8, truth.basisData.getChannel(0), recon.basisData.getChannel(0));
  std::cout << "depth:" << stats  << std::endl;
  std::cout << "\n";


  // Output the stats for every frame.
  image_stats_t flowStatsAverage;
  image_stats_t posStatsAverage;

  if (zeroDisplacements)
    recon.basisData = recon.basisData * 0;

  for (int i=0; i<numTimes; i++) {
    truth.setTime(i);
    recon.setTime(i);

    nacb::Imagef DtruthWorld, DreconWorld;
    nacb::Imagef Dtruth = truth.basis->getDisplacementMap(truth.basisData, i);
    nacb::Imagef Drecon = recon.basis->getDisplacementMap(recon.basisData, i);
    nacb::Imagef posTruth = convertDisplacementToPosition(truth.getBaseGeometry(), Dtruth, DtruthWorld);
    nacb::Imagef posRecon = convertDisplacementToPosition(recon.getBaseGeometry(), Drecon, DreconWorld);

    // For reference time compute the depths
    if (i == 0) {
      std::cout << "Projected depth.\n";

      depthTruth = convertPositionToDisplacement(truth.getBaseGeometry(), posTruth);
      depthRecon = convertPositionToDisplacement(recon.getBaseGeometry(), posRecon);

      double n, x;
      depthTruth.getChannel(0).getRange(n, x);
      std::cout << "Truth depth range:" << n << " " << x << std::endl;
      
      depthRecon.getChannel(0).getRange(n, x);
      std::cout << "Recon depth range:" << n << " " << x << std::endl;
      
      depthTruth.getChannel(0).getNormalizedImage().write("/tmp/t2.png");
      depthRecon.getChannel(0).getNormalizedImage().write("/tmp/r2.png");
      
      image_stats_t truth_zero = disp_compare(mask8, depthTruth, depthRecon);
      std::cout << "truth (zero):" << truth_zero << std::endl;
      
      image_stats_t recon_zero = disp_compare(mask8, depthRecon, depthTruth);
      std::cout << "recon (zero):" << recon_zero << std::endl;
      
      image_stats_t stats = disp_compare(mask8, depthTruth, depthRecon);
      std::cout << "depth:" << stats  << std::endl;
      std::cout << "\n";
    }

    image_stats_t flowStats = flow_compare(mask8, DtruthWorld, DreconWorld);

    std::cout << "flow" << i << ":" << flowStats << std::endl;

    // Now displace the points
    image_stats_t posStats = flow_compare(mask8, posTruth, posRecon);

    std::cout << "posn" << i << ":" << posStats << std::endl;
    std::cout << "\n";

    if (i == 0) {
      flowStatsAverage = flowStats;
      posStatsAverage = posStats;
    }
    else {
      flowStatsAverage += flowStats;
      posStatsAverage += posStats;
    }
  }

  flowStatsAverage /= numTimes;
  posStatsAverage /= numTimes;

  std::cout << "averages." << std::endl;
  std::cout << "flow:" << flowStatsAverage << std::endl;
  std::cout << "posn:" << posStatsAverage << std::endl;


  std::cout << "Loading seqs:" << av[optind] << std::endl;

  std::vector<DisparitySequence> seqs = 
    DisparitySequence::loadSequences(av + optind, ac - optind);

  if (!seqs.size()) exit(1);

  imageComparison(seqs, recon, truth, numTimes);

  return 0;
}
