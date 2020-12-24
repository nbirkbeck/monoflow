#include "DisparityMap.h"
#include "DisparitySequence.h"
#include "discrete.h"
#include "utils.h"

#include <nmisc/commandline.h>

#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <vector>


namespace bfs=boost::filesystem;


int main(int ac, char * av[]){
  nacb::CommandLine cline;
  int downsample = 1;
  bool solveDiscreteFlow = false;
  double dmin = 5.0, dmax = 9.0;
  double smoothWeight = kSmoothWeight;
  int numDisparityIterations = 2;
  int numDisparityLevels = 60;
  int nclosest = -1;

  std::string dispName = "";
  std::string oname = "/tmp/discrete-%04d.rfi";
  std::string odepth = "";
  std::string debugRoot = "";

  cline.registerOption("downsample", "Number of times to perform downsampling.", &downsample, 'd');
  cline.registerOption("dmin", "Min disparity", &dmin, 0);
  cline.registerOption("dmax", "Min disparity", &dmax, 0);

  cline.registerOption("disp", "Input disparity maps.", &dispName, 0);

  cline.registerOption("ndisp", "Number of disparity levels.", &numDisparityLevels, 0);
  cline.registerOption("diters", "Number of disparity iterations.", &numDisparityIterations, 0);
  cline.registerOption("odisp", "Output name for disparity maps.", &oname, 0);
  cline.registerOption("odepth", "Output name for depth maps.", &odepth, 0);
  cline.registerOption("smooth", "Smoothness weight.", &smoothWeight, 0);
  cline.registerOption("nclosest", "Use the n-closest (in time) images for computation", &nclosest, 0);
  
  cline.registerOption("root", "Debug root", &debugRoot);
  cline.parse(ac, av);
  
    
  if (debugRoot.size())
    bfs::create_directory(debugRoot);

  std::vector<SeqType> seqs = SeqType::downsample(SeqType::loadSequences(av + optind, ac - optind), downsample);
  std::vector<DiscreteDisparityMapInt> disps;

  if(nclosest <= 0 || nclosest >= (int)seqs.size())
    nclosest = seqs.size();

  // Load a saved disparity; allows for re-running at higher resolutions.
  if(dispName.size()){

    printf("Loading disparity from %s\n", dispName.c_str());

    for(int i=0; i<(int)seqs.size(); i++){
      nacb::Imagef d((boost::format(dispName) % i).str().c_str());
      
      d = d.resize(seqs[0].image.w, seqs[0].image.h);
      
      DiscreteDisparityMapInt dmap(d.w, d.h, 1.0/dmax, 1.0/dmin, numDisparityLevels);
      dmap.setFromDisparity(d);

      disps.push_back(dmap);
    }
  }

  if(numDisparityIterations){
    for(int its=0; its<numDisparityIterations; its++){
      printf("Solving iteration: %d\n", its);

      solve_disparity(seqs, 1.0/dmax, 1.0/dmin, numDisparityLevels, disps, smoothWeight, nclosest, downsample);
    
      /**
	 Only output if the debug root is not-null.
      */
      if (debugRoot.size()) {
	std::string dir = (boost::format(debugRoot + "/ndi-%02d") % its).str();
	bfs::create_directory(dir);

	write_all_warped_images(seqs, disps, (dir + "/w-%02d-%02d.png").c_str());
	
	for(int j=0; j<(int)seqs.size(); j++){
	  int w = seqs[j].image.w;
	  int h = seqs[j].image.h;
	  std::string fname = (boost::format("depth-%d.rfi") % j).str();
	  disparity_to_depth(disps[j].getMappedDisparity()).resize(w<<downsample, h<<downsample).write((dir + "/" + fname).c_str());
	}
      }

      // Compare disparities.
      for(int i=0; i<(int)disps.size(); i++){
	if(seqs[i].disp.nchannels > 0)
	  disp_compare(seqs[i].image, seqs[i].disp, disps[i].getMappedDisparity());
      }
    }

    /* 
       The discrete flow was a pretty bad idea.  Too many variables are needed to
       describe the motion.
    */
    if(solveDiscreteFlow){
      for(int its=0; its<1; its++){
	std::vector<DiscreteFlowMap> flowMaps = solve_flow(seqs, disps, 0.05, 8, 8, 0);
      
	write_all_warped_images(seqs, disps, flowMaps, "/tmp/flow-%02d-%02d.png");
      }
    }
  }

  // Write the output files (if required)
  if (oname.size() || odepth.size()) {
    for(int i=0; i<(int)disps.size(); i++) {
      int w = seqs[i].image.w;
      int h = seqs[i].image.h;
      nacb::Imagef disp = disps[i].getMappedDisparity();
      disp = disp.resize(w<<downsample, h<<downsample);
      
      if(oname.size())
	disp.write((boost::format(oname) % i).str().c_str());

      if(odepth.size())
	disparity_to_depth(disp).write((boost::format(odepth) % i).str().c_str());
    }
  }

  return 0;
}
