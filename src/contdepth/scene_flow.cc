/**
   This is meant as a test of scene flow

   Might be useful to have separate constraints on the depth.
*/
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <nmath/vec3.h>
#include <nmisc/commandline.h>


#include <vector>
#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "autorecon/stereo/conversions.h"
#include "autorecon/stereo/flow.h"
#include "EpipolarLine.h"
#include "DisparityMap.h"
#include "FlowMap.h"

#include "variational.h"
#include "variational-sf.h"
#include "discrete.h"
#include "DisparitySequence.h"


void scene_flow_time(std::vector<SeqType> & seqs, 
		     std::vector<SeqType> & seqs_next, 
		     std::vector<nacb::Imagef> & flows,
		     double imageCostWeight,
		     int downsample){
  double h = 1;
  
  for(int source = 0; source<seqs.size(); source++){
    int upw = seqs[source].image.w * (1<<downsample);
    int uph = seqs[source].image.h * (1<<downsample);
    
    // Test out the variational approximation
    std::vector<DisparityOffset> sfwarps;
    std::vector<DisparityOffset> sfwarp_invs;

    std::vector<nacb::Imagef> images;

    // The first warps do not have offsets (hence the 0) below
    images.push_back(seqs_next[source].image);
    sfwarps.push_back(DisparityOffset(seqs[source].A, seqs[source].E, 
				      seqs[source].A, seqs[source].E, 1));
    sfwarp_invs.push_back(DisparityOffset(seqs[source].A, seqs[source].E, 
					  seqs[source].A, seqs[source].E, -1));
      

    for(int i=0; i<seqs.size(); i++){
      if(i == source)
	continue;

      // The first warps do not have offsets (hence the 0) below
      images.push_back(seqs[i].image);
      sfwarps.push_back(DisparityOffset(seqs[source].A, seqs[source].E, 
					seqs[i].A, seqs[i].E, 0));

      sfwarp_invs.push_back(DisparityOffset(seqs[i].A, seqs[i].E, 
					    seqs[source].A, seqs[source].E, 0));
       
      images.push_back(seqs_next[i].image); //im1 = warp_image(seqs, disps, 0, 1);

      // Assume the sequences are ordered
      sfwarps.push_back(DisparityOffset(seqs[source].A, seqs[source].E, 
					seqs[i].A, seqs[i].E, 1));

      sfwarp_invs.push_back(DisparityOffset(seqs[i].A, seqs[i].E, 
					    seqs[source].A, seqs[source].E, -1));
    }
    VariationalSceneFlowProblem flow(seqs[source].image, sfwarps, 
				     images, flows[source], imageCostWeight,
				     VariationalSceneFlowProblem::defaultBeta, h, h);

    nacb::Imagef residuals;
    mg::wcycle_solver solver(3, 2, 7, 7, true);

    flow.iteration(residuals, false).print();
    solver.solve(flow);
    flow.iteration(residuals, false).print();

    nacb::Imagef soln = flow.getCombinedResult();
    flows[source] = soln;

    
    (flows[source].getChannel(0)).getNormalizedImage().write("/tmp/flow-d.png");
    flows[source].getChannel(1).getNormalizedImage().write("/tmp/flow-u.png");
    flows[source].getChannel(2).getNormalizedImage().write("/tmp/flow-v.png");
    flows[source].getChannel(3).getNormalizedImage().write("/tmp/flow-w.png");
    
    std::string dirname = "/tmp/";
    // Warp all the images
    for(int j=0; j<=images.size(); j++){
      std::string fname = dirname + (boost::format("sf-%02d-%02d.png") % source % j).str();
	  
      if(j == 0){
	seqs[source].image.write(fname.c_str());
	continue;
      }  
      
      warp_image(seqs[source].image.w, seqs[source].image.h, images, 
		 sfwarps, flows[source], j-1).write(fname.c_str());


    }
    std::string fname = dirname + (boost::format("uv-flow-%d.flo") % source).str();
    flow_save_flo(project_flow(seqs[source].A, seqs[source].E, flows[source]), fname.c_str());
  }  
}


/**
   
 */
int main(int ac, char *av[]){
  nacb::CommandLine cline;
  int downsample = 1;
  int variationalDisparity = 1;
  int variationalFlow = 1;
  int useTruth = 0;
  int vdisp_max = 0;
  double flowImageCostWeight = 0.15;
    
  std::string depthName = "", dispName = "", flowName = "";
  std::string debugRoot = "/tmp/cout";
  std::string odisp = "", odepth = "", oflow = "";

  int index1 = 0;
  int index2 = 20;

  cline.registerOption("downsample", "Number of times to perform downsampling.", &downsample, 'd');
  cline.registerOption("vdisp_dsmax", "Maximum downsample to reach for disparity.", &vdisp_max, 0);
  cline.registerOption("vdisp", "Do variational disparity", &variationalDisparity, 0);
  cline.registerOption("vflow", "Do variational flow", &variationalFlow, 0);
  cline.registerOption("use_truth", "Use truth", &useTruth, 0);

  cline.registerOption("flow", "Flow input name", &flowName, 0);
  cline.registerOption("disp", "Name for isparity maps.", &dispName, 0);
  cline.registerOption("depth", "Name for input depth maps (instead of disparity).", &depthName, 0);

  cline.registerOption("index1", "The index of the first image to load.", &index1, 0);
  cline.registerOption("index2", "The index of the first image to load.", &index2, 0);

  cline.registerOption("root", "Debug root", &debugRoot, 0);
  cline.registerOption("odisp", "Output name for disparity maps.", &odisp, 0);
  cline.registerOption("odepth", "Output name for depth maps.", &odepth, 0);
  cline.registerOption("oflow", "Output name for flows.", &oflow, 0);
  cline.registerOption("fimg", "Image weight", &flowImageCostWeight, 0);
  cline.parse(ac, av);

  std::vector<SeqType> seqs_full_res = DisparitySequence::loadSequences(av + optind, ac - optind);
  std::vector<SeqType> seqs, seqs_next;
  std::vector<SeqType> seqs_next_full_res = DisparitySequence::loadSequences(av + optind, ac - optind);
  
  if(downsample){
    seqs = DisparitySequence::downsample(seqs_full_res, downsample);
    seqs_next = DisparitySequence::downsample(seqs_next_full_res, downsample);
  }
  else {
    seqs = seqs_full_res;
    seqs_next = seqs_next_full_res;
  }
  printf("Loading.\n");

  for(int i=0; i<seqs.size(); i++){
    int w = seqs[i].image.w;
    int h = seqs[i].image.h;

    seqs_full_res[i].load((int)index1);
    seqs[i].load((int)index1);

    seqs_next_full_res[i].load(index2);
    seqs_next[i].load(index2);
    printf("after\n");

    // Move this crap
    if(downsample){
      seqs[i].image = seqs[i].image.resize(w, h);
      seqs_next[i].image = seqs_next[i].image.resize(w, h);
    }
  }

  double imageCostWeight = flowImageCostWeight; //0.05;

  std::vector<nacb::Imagef> flows;
  std::vector<nacb::Imagef> disps;

  if(dispName.size()){    
    printf("Using disps\n");

    for(int i=0; i<seqs.size(); i++) {
      disps.push_back(nacb::Imagef((boost::format(dispName) % i).str().c_str()));
      disps.back() = disps.back().resize(seqs[i].image.w, seqs[i].image.h);
    }
  }
  else if(depthName.size()){
    // Load in the depth maps 
    printf("Using depth maps.\n");

    for(int i=0; i<seqs.size(); i++) {
      Imagef depth = nacb::Imagef((boost::format(depthName) % i).str().c_str());
      printf("Loaded depth:%d  %dx%d (ignoring < 0)\n", i, depth.w, depth.h);
      depth = (depth < 0) + (depth > 0)*depth;
      disps.push_back(disparity_to_depth(depth));
      disps.back() = disps.back().resize(seqs[i].image.w, seqs[i].image.h);
    }
  }
  else
    assert(0);

  for(int i=0; i<disps.size(); i++){
    nacb::Imagef flow(disps[i].w, disps[i].h, 4);
    flow = 0;
    flow.setChannel(0, disps[i]);
    flows.push_back(flow);
  }

  printf("resizing images.\n");
  disps = resize_images(disps, seqs[0].image.w, seqs[0].image.h);
  flows = resize_images(flows, seqs[0].image.w, seqs[0].image.h);
  printf("done resize\n");

  while(downsample >= vdisp_max){

    for(int its=0; its<3; its++)
      scene_flow_time(seqs, seqs_next, flows, imageCostWeight, downsample);

    downsample--;

    if(downsample > 0){
      seqs = DisparitySequence::downsample(seqs_full_res, downsample);
      seqs_next = DisparitySequence::downsample(seqs_next_full_res, downsample);
    }
    else if(downsample == 0){
      seqs = seqs_full_res;
      seqs_next = seqs_next_full_res;
    }

    flows = resize_images(flows, seqs[0].image.w, seqs[0].image.h);
  }

  save_depth_and_flow(flows, odepth, oflow);

  return 0;
}
