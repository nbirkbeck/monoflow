/**
   A function to compare the results of the disparity
   estimation (including warped images and to ground truth).
 */
#include "DisparitySequence.h"
#include "EpipolarLine.h"
#include "utils.h"
#include <nmisc/commandline.h>
#include <nimage/image.h>
#include <vector>


typedef double ImageStats;

std::vector<image_stats_t> warp_image_compare(const std::vector<SeqType> & seqs,
					   const nacb::Imagef & disp,
					   int source,
					   int nclosest = -1){
  std::vector<image_stats_t> imageStats;

  int low, high;
  get_closest_range(nclosest, seqs.size(), source, low, high);

  for(int i=low; i<high; i++){
    EpipolarLine eline(seqs[source].A, seqs[source].E,
		       seqs[i].A, seqs[i].E);
    nacb::Image8 mask = seqs[source].image.getChannel(3)>100;

    if(i != source){
      nacb::Image8 image2 = warp_image(seqs[i].image, eline, disp);
      imageStats.push_back(image_compare(mask, seqs[source].image, image2));
    }
  }
  return imageStats;
}


std::vector<image_stats_t> warp_image_compare(const std::vector<SeqType> & seqs,
					   const nacb::Imagef & disp,
					   const nacb::Imagef & flow,
					   int source,
					   int nclosest = -1){
  std::vector<image_stats_t> imageStats;

  int low, high;
  get_closest_range(nclosest, seqs.size(), source, low, high);
  
  for(int i=low; i<high; i++){
    DisparityOffset warp(seqs[source].A, seqs[source].E,
			 seqs[i].A, seqs[i].E, double(i - source));
    nacb::Image8 mask = seqs[source].image.getChannel(3)>100;

    if(i != source){
      nacb::Image8 image2 = warp_image(seqs[i].image, warp, disp, flow);
      imageStats.push_back(image_compare(mask, seqs[source].image, image2));
    }
  }
  return imageStats;
}



std::ostream& operator<<(std::ostream & o, const std::vector<image_stats_t> & stats){
  ImageStats mean = 0;
  
  for(int i=0; i<(int)stats.size(); i++)
    mean += stats[i].mean;

  mean /= stats.size();

  o << mean << " ";
  
  o << "[";
  
  for(int i=0; i<(int)stats.size(); i++)
    o << stats[i].mean << " ";

  o << "]";
  return o;
}


int main(int ac, char * av[]){
  nacb::CommandLine cline;
  std::string dispName = "", depthName = "";
  std::string flowName = "";
  std::string flowWarpName = "";
  std::string dispWarpName = "";

  int nclosest = -1;

  cline.registerOption("disp", "Name for disparity maps.", &dispName, 0);
  cline.registerOption("flow", "Name for flows.", &flowName, 0);
  cline.registerOption("depth", "Name for depth maps.", &depthName, 0);
  cline.registerOption("nclosest", "Number of closest images to use.", &nclosest, 0);
  cline.registerOption("disp-warp", "Name of images for disp-warp (should contain two %d's)", &dispWarpName, 0);
  cline.registerOption("flow-warp", "Name of the images for flow-warp (should contain two %d's)", &flowWarpName, 0);
  cline.parse(ac, av);
  
  std::vector<SeqType> seqs = DisparitySequence::loadSequences(av + optind, ac - optind);
  std::vector<nacb::Imagef>  disps, flows;
  

  printf("%s\n\n\n", dispName.c_str());

  for(int i=0; i<(int)seqs.size(); i++){    
    if(dispName.size())
      disps.push_back((boost::format(dispName) % i).str().c_str());
    
    else if(depthName.size())
      disps.push_back(disparity_to_depth((boost::format(depthName) % i).str().c_str()));

    if(flowName.size())
      flows.push_back((boost::format(flowName) % i).str().c_str());
  }

  for(int i=0; i<(int)seqs.size(); i++){
    printf("sequence:%d\n", i);
    std::vector<image_stats_t> stats;

    if(flows.size())
      stats = warp_image_compare(seqs, disps[i], flows[i], i, nclosest);
    else 
      stats = warp_image_compare(seqs, disps[i], i, nclosest);

    if(seqs[i].disp.w >= 1){
      image_stats_t s = disp_compare(seqs[i].image, seqs[i].disp, disps[i]);
      
      std::cout << "  disp:" << s << "\n";

      image_stats_t sd = disp_compare(seqs[i].image, disparity_to_depth(seqs[i].disp), disparity_to_depth(disps[i]));
      
      std::cout << "  depth:" << sd << "\n";
    }

    if(flows.size() && seqs[i].flow.w >= 1)
      std::cout << "  disp:" << flow_compare(seqs[i].image, seqs[i].flow, flows[i]) << "\n";;

    std::cout << "  image(rms): " << stats << "\n"; 
  }

  if (dispWarpName.size() && disps.size()) {
    write_all_warped_images(seqs, disps, dispWarpName.c_str(), nclosest, WARP_ALPHA_SOURCE);
  }

  if (flowWarpName.size() && flows.size()) {
    write_all_warped_images(seqs, disps, flows, flowWarpName.c_str(), nclosest, WARP_ALPHA_SOURCE);
  }

  return 0;
}
