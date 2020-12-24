/**
  TODO:
   Subpixel refinement?
   
   How to sample the cost.
   Analytic derivative for displacement consistency (these are okay now),
     -There are still some nans (sometimes).
  
  Variational disparity: 
   -Allow for a subset of the images to be used (e.g., best matching to start with).   
     -Filter this spatially maybe?  
     -Maybe use the disparity agreement for this part too.
     
  Flow:
   -Add consistency term... it is pretty damn complicated, me thinks.
   -
  -Functions to compute the cost terms (and display them, compare them to ground truth).
  How to initialize? Coarse to fine, discrete BP, some combination?

  Smoothness parameters!
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

#include "EpipolarLine.h"
#include "DisparityMap.h"
#include "FlowMap.h"

#include "variational.h"
#include "variational-sf.h"
#include "discrete.h"
#include "DisparitySequence.h"


using namespace nacb;
namespace bfs=boost::filesystem;

const double kDisparityDispConstWeightGood = 0.0001;
const double kDispImageCostWeight = 0.25;

template <bool UseDisparity>
void variational_disparity_output(const std::string & dirname_,
				  const VariationalDispProblem<UseDisparity> & vari,
				  const std::vector<SeqType> & seqs,
				  const std::vector<EpipolarLineTemplate<UseDisparity> > & elines,
				  const std::vector<nacb::Imagef> & images,
				  const std::vector<nacb::Imagef> & disps,
				  int source,
				  int downsample){
  int upw = seqs[source].image.w * (1<<downsample);
  int uph = seqs[source].image.h * (1<<downsample);

  std::string dirname = dirname_ + "/";
  const nacb::Imagef & soln = disps[source];

  // Print out some stats regarding the update estimation.
  bool printRange = !false;
  
  double n, x;
  soln.getRange(n, x);
  
  if(printRange)
    printf("range: %lf %lf\n", n, x);
    
  std::string ofile = dirname + "vd-%02d";
  std::string fname = (boost::format(ofile) % source).str();
  ((soln - n)*(1.0/(x-n))).write((fname + ".png").c_str());
  
  if (UseDisparity)
    disparity_to_depth(soln).resize(upw, uph).write((fname + ".rfi").c_str());
  else 
    soln.resize(upw, uph).write((fname + ".rfi").c_str());

  // FIXME: Didn't used to get the range here...
  vari.D_0.getRange(n, x);
  if (printRange)
    printf("d0 range: %lf %lf\n", n, x);

  ((vari.D_0.resize(upw, uph) - n)*(1.0/(x-n))).write((fname + "_init.png").c_str());
  
  if(seqs[source].disp.nchannels) {
    if (UseDisparity)
      disp_compare(seqs[source].image, seqs[source].disp, soln);
    else 
      fprintf(stderr, "Warning: not comparing to the ground truth disparity (as we are using depth representation.\n");
  }
}

template <bool UseDisparity> 
struct vdisp_output {
  typedef  boost::function< void (const VariationalDispProblem<UseDisparity> & vari,
				  const std::vector<EpipolarLineTemplate<UseDisparity> > & elines,
				  const std::vector<nacb::Imagef> & images,
				  const std::vector<Imagef> & disps,
				  int source,
				  int downsample)> func_t;
};



template <bool UseDisparity>
void variational_disparity(std::vector<SeqType> & seqs,
			   std::vector<Imagef> & disps,
			   int downsample,
			   double imageCostWeight = kDispImageCostWeight,
			   double dispConstWeight = kDisparityDispConstWeightGood,
			   const typename vdisp_output<UseDisparity>::func_t & output_func = 
			         typename vdisp_output<UseDisparity>::func_t(),
			   int nclosest = -1){
  double smoothnessWeight = 10;
  
  int nuse = std::min(2*nclosest, (int)seqs.size());  
  imageCostWeight = imageCostWeight/nuse;

  /*
  for(int i=0; i<seqs.size(); i++)
    disps[i] = disps[i].boxFilter(10);
  */

  // testing the plane.
  /*
  for(int i=0; i<seqs.size(); i++){
    disps[i] = 0.12;
  }
  */
  
  vector<vector<Vec4d> > results;
  for(int i=0; i<(int)seqs.size(); i++){
    results.push_back(std::vector<Vec4d>());
  }

  //int numits = 3;
  //for(int its=0; its<numits; its++)
  {  
    for(int source=0; source<(int)seqs.size(); source++){
      std::vector<EpipolarLineTemplate<UseDisparity> > elines, eline_invs;
      std::vector<nacb::Imagef> images, disps_few;
    
      nacb::Imagef dsource = disps[source];
      nacb::Imagef weight = weight_from_alpha(seqs[source].image);

      int lo, hi;
      get_closest_range(nclosest, seqs.size(), source, lo, hi);

      for(int i=lo; i<hi; i++){
	if(i == source)
	  continue;
      
	images.push_back(seqs[i].image); //im1 = warp_image(seqs, disps, 0, 1);
	elines.push_back(EpipolarLineTemplate<UseDisparity>(seqs[source].A, seqs[source].E, seqs[i].A, seqs[i].E));

	eline_invs.push_back(EpipolarLineTemplate<UseDisparity>(seqs[i].A, seqs[i].E, seqs[source].A, seqs[source].E));
	disps_few.push_back(disps[i]);
      }
    
      VariationalDispProblem<UseDisparity> vari(seqs[source].image, elines, 
						images, dsource, imageCostWeight, smoothnessWeight);

      
      {
	double n, x;
	weight.getRange(n, x);
	((weight - n)*(1.0/(x-n))).save("/tmp/weight.png");
	printf("weight: %f %f\n", n, x);
      }

      // This doesn't work all that well.
      if(fabs(dispConstWeight) > 0){
	vari.addDisplacementConsistency(elines, eline_invs, disps_few, dispConstWeight);
      }

      // Evaluate actual energy
      double S_energy = VariationalDispProblem<UseDisparity>::smoothnessEnergy(dsource);
      double G_energy = VariationalDispProblem<UseDisparity>::geometryDataEnergy(dsource, weight, elines, eline_invs, disps_few);
      double I_energy = VariationalDispProblem<UseDisparity>::imageDataEnergy(dsource, seqs[source].image, weight, elines, images);
      
      results[source].push_back(Vec4d(G_energy*dispConstWeight + S_energy*smoothnessWeight + I_energy*imageCostWeight,
				      S_energy, G_energy, I_energy));
      
      printf("Total energy(%d): %lf (%f %f %f)\n", 
	     source, (dispConstWeight>0) ? G_energy*dispConstWeight : 0.0
	     + S_energy*smoothnessWeight + I_energy*imageCostWeight,
	     G_energy, S_energy, I_energy);

      nacb::Imagef residuals;
      mg::wcycle_solver solver(3, 3, 10, 10, true);

    
      // Some presmoothing to debug some nans....
      for(int i=0; i<0; i++){
	vari.iteration(residuals, true).print();
      }

      vari.iteration(residuals, false).print();
      solver.solve(vari);
      vari.iteration(residuals, false).print();    
      
    
      nacb::Imagef soln = vari.getCombinedResult();
      disps[source] = soln;
      dsource = soln;

      output_func(vari, elines, images, disps, source, downsample);

      if(fabs(dispConstWeight) > 0){
	VariationalDispProblem<UseDisparity> vari2(seqs[source].image, elines, 
						   images, dsource, imageCostWeight, smoothnessWeight);
      
	vari2.addDisplacementConsistency(elines, eline_invs, disps_few, dispConstWeight);
      }

      G_energy = VariationalDispProblem<UseDisparity>::geometryDataEnergy(dsource, weight, elines, eline_invs, disps_few);
      S_energy = VariationalDispProblem<UseDisparity>::smoothnessEnergy(dsource);
      I_energy = VariationalDispProblem<UseDisparity>::imageDataEnergy(dsource, seqs[source].image, weight, elines, images);

      results[source].push_back(Vec4d(G_energy*dispConstWeight + 
				      S_energy*smoothnessWeight + 
				      I_energy*imageCostWeight,
				      S_energy, G_energy, I_energy));
      G_energy = VariationalDispProblem<UseDisparity>::geometryDataEnergy(dsource, weight, elines, eline_invs, disps_few);
      printf("Total energy (af): %d (%f %f %f %f)\n", 
	     source, (dispConstWeight>0 ? G_energy*dispConstWeight : 0)
	     + S_energy*smoothnessWeight + I_energy*imageCostWeight,
	     G_energy, S_energy, I_energy);
    }
  }

  for(int i=0; i<(int)results.size(); i++){
    for(int j=0; j<(int)results[i].size(); j++){
      std::cout << results[i][j] << "\n";
    }
    std::cout << "\n\n\n";
  }
}

template <class DisparityOffsetType>
void variational_scene_flow_output(const std::string & dirname_,
				   const VariationalSceneFlowProblem<DisparityOffsetType::UseDisparity> & problem,
				   std::vector<SeqType> & seqs,
				   const std::vector<DisparityOffsetType> & sfwarps,
				   std::vector<nacb::Imagef> & images,
				   std::vector<Imagef> & flows,
				   int source,
				   int downsample){
  bool UseDisparity = DisparityOffsetType::UseDisparity;
  std::cerr << "scene_flow_output to " << dirname_ <<  "( UseDisparity:" << UseDisparity << ")\n";

  int upw = seqs[source].image.w << downsample;
  int uph = seqs[source].image.h << downsample;
  
  nacb::Imagef soln = flows[source];
  std::string num = (boost::format("-%d") % source).str();
  std::string imnames[4]= {"dd", "du", "dv", "dw"};  
  std::string dirname = dirname_ + "/";

  for(int j=0; j<4; j++){		      
    std::cout << "Writing channel: " << j << "\n";
    
    double n, x;
    nacb::Imagef im = soln.getChannel(j);
    
    im.getRange(n, x);	
    im.getNormalizedImage().write((dirname + imnames[j] + num + ".png").c_str());
    
    if(j >= 1) {
      nacb::Imagef chan;
      chan = seqs[source].flow.getChannel(j - 1);
      ((chan - n)*(1.0/(x - n))).write((dirname + imnames[j] + num + "_truth.png").c_str());
    }
    
    printf("%d range: %f %f\n", j, n, x);
  }
  
  nacb::Imagef flowOnly(soln.w, soln.h, 3);
  flowOnly.setChannel(0, soln.getChannel(1));
  flowOnly.setChannel(1, soln.getChannel(2));
  flowOnly.setChannel(2, soln.getChannel(3));
  flowOnly.resize(upw, uph).write((dirname + std::string("flow") + num + ".rfi").c_str());
  
  std::string d_flow = dirname + "d-flow" + num + ".rfi";
  if (UseDisparity) {
    disparity_to_depth(soln.getChannel(0)).resize(upw, uph).write(d_flow.c_str());
  }
  else {
    soln.getChannel(0).resize(upw, uph).write(d_flow.c_str());
  }
  
  // Warp all the images
  for(int j=0; j<=(int)images.size(); j++){
    std::string fname = dirname + (boost::format("sf-%02d-%02d.png") % source % j).str();
    
    if(j == 0){
      seqs[source].image.write(fname.c_str());
      continue;
    }
    warp_image(seqs[source].image.w, seqs[source].image.h, images, 
	       sfwarps, soln, j-1).write(fname.c_str());
  }
}


const double kSceneFlowDispContWeightGood = 0.0005; // A good number for the disp cont weight.
const double kSceneFlowImageCostWeight = 0.1;

template <bool UseDisparity> 
class scene_flow_output {
public:
  typedef boost::function< void (const VariationalSceneFlowProblem<UseDisparity> &,
				 const std::vector<DisparityOffsetBase<UseDisparity> > & sfwarps,
				 std::vector<nacb::Imagef> & images,
				 std::vector<Imagef> & flows,
				 int source,
				 int downsample)> func_t;
};



template <class DispOffsetType>
void variational_scene_flow(std::vector<SeqType> & seqs,
			    std::vector<Imagef> & flows,
			    int downsample = 0,
			    double imageCostWeight = kSceneFlowImageCostWeight,
			    double flowBeta = VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::defaultBeta.y,
			    double dispContWeight = 0,
			    const typename scene_flow_output<DispOffsetType::UseDisparity>::func_t & output_func = 
			          typename scene_flow_output<DispOffsetType::UseDisparity>::func_t(),
			    int nclosest = -1){
  //const bool testDispCont = false;
  const bool constantDisparity = false; //Used for debugging the flow only (does it decrease score).
  double h = 1; // You would think that this should take into account the downsampling, so regulariz is consistent, but that doesn't work.

  // Normalize the cost by the number of sequences.
  imageCostWeight /= seqs.size();

  Vec2d beta = VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::defaultBeta;
  beta.y = flowBeta;

  for(int source=0; source<(int)seqs.size(); source++){
    // Test out the variational approximation
    std::vector<DispOffsetType> sfwarps;
    std::vector<DispOffsetType> sfwarp_invs;
    std::vector<nacb::Imagef> images;
    std::vector<nacb::Imagef> neigh_flows;

    int lo, hi;
    get_closest_range(nclosest, seqs.size(), source, lo, hi);
    
    for(int i=lo; i<hi; i++){
      if(i == source)
	continue;

      images.push_back(seqs[i].image); //im1 = warp_image(seqs, disps, 0, 1);

      // Assume the sequences are ordered
      sfwarps.push_back(DispOffsetType(seqs[source].A, seqs[source].E, 
				       seqs[i].A, seqs[i].E, double(i - source)));
      
      sfwarp_invs.push_back(DispOffsetType(seqs[i].A, seqs[i].E, 
					   seqs[source].A, seqs[source].E, double(source - i)));

      neigh_flows.push_back(flows[i]);
    }
    VariationalSceneFlowProblem<DispOffsetType::UseDisparity>
      flow(seqs[source].image, sfwarps, 
	   images, flows[source], imageCostWeight,
	   beta, h, h);

    if(dispContWeight > 0){
      flow.addDisplacementConsistency(sfwarps, sfwarp_invs, neigh_flows, dispContWeight);
    }

    if(constantDisparity)
      flow.addDisparityConstraint();

    // Get initial residuals
    {
      double before = 
	VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::imageDataEnergy(flows[source], seqs[source].image, 
										   flow.weight, sfwarps, images);

      double geomBefore =
	VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::geometricDataEnergy2D(flows[source], seqs[source].image, 
											 flow.weight, sfwarps, sfwarp_invs,
											 images, neigh_flows);
      double geom3D = 
	VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::geometricDataEnergy3D(flows[source], seqs[source].image, 
											 flow.weight, sfwarps, sfwarp_invs,
											 images, neigh_flows);
      
      std::cout << "Before: " << before << " " << geomBefore << " " << geom3D << "\n";
    }

    nacb::Imagef residuals;
    mg::wcycle_solver solver(3, 2, 7, 7, true);

    // Just run a couple iterations to see the residuals.
    for(int i=0; i<3; i++){
      flow.iteration(residuals, true).print();
    }
    
    flow.iteration(residuals, false).print();
    solver.solve(flow);
    flow.iteration(residuals, false).print();

    nacb::Imagef soln = flow.getCombinedResult();
    flows[source] = soln;


    output_func(flow, sfwarps, images, flows, source, downsample);
    
    double after = 
      VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::imageDataEnergy(flows[source], seqs[source].image, 
										 flow.weight, sfwarps, images);
    double geomAfter =
      VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::geometricDataEnergy2D(flows[source], seqs[source].image, 
										       flow.weight, sfwarps, sfwarp_invs,
										       images, neigh_flows);
    double geom3D = 
      VariationalSceneFlowProblem<DispOffsetType::UseDisparity>::geometricDataEnergy3D(flows[source], seqs[source].image, 
										       flow.weight, sfwarps, sfwarp_invs,
										       images, neigh_flows);
    
    std::cout << "After: " << after << " " << geomAfter  << " " << geom3D << "\n\n";
  }
}


void write_disparities_and_depths(const std::vector<Imagef> & disps, 
				  const std::string & odisp,
				  const std::string & odepth,
				  int w, int h){
  for(int i=0; i<(int)disps.size(); i++){
    Imagef  disp = disps[i].getChannel(0).resize(w, h);
    
    if(odisp.size())
      disp.write((boost::format(odisp) % i).str().c_str());
    
    if(odepth.size())
      disparity_to_depth(disp).write((boost::format(odepth) % i).str().c_str());
  }
}

template <bool UseDisparity>
void variational_disparity_main(std::vector<SeqType>& seqs_full_res,
				std::vector<SeqType>& seqs,
				std::vector<Imagef>& disps,
				double dispImageCostWeight,
				double dispConstWeight,
				int variationalDisparity,
				int nclosest,
				int downsample,
				int vdisp_max,
				const std::string& odisp,
				const std::string& odepth,
				const std::string& debugRoot)
{
  std::cout << "-------------------------------------\n\n";
  int downsample_local = downsample;

  if (!UseDisparity) {
    for (int i = 0; i < (int)disps.size(); ++i) {
      disps[i] = disparity_to_depth(disps[i]);
    }
  }

  while(downsample_local >= vdisp_max){
    typename vdisp_output<UseDisparity>::func_t output_func;
    
    if (debugRoot.size()) {
      std::string output_dir = "";
      
      output_dir = (boost::format(debugRoot + "/vdisp/r%d") % downsample_local).str();
      bfs::create_directory(debugRoot + "/vdisp");
      bfs::create_directory(output_dir);
      
      output_func = boost::bind(::variational_disparity_output<UseDisparity>, output_dir, _1, 
				seqs, _2, _3, _4, _5, _6);
    }
    
    printf("solving disparity for downsample=%d\n", downsample_local);
    
    for(int i=0; i<variationalDisparity; i++)
      variational_disparity<UseDisparity>(seqs, disps, downsample_local, dispImageCostWeight * sqrt(double(1 << downsample_local)), dispConstWeight, output_func, nclosest);
    
    downsample_local--;
    
    if(downsample_local > 0){
      seqs = SeqType::downsample(seqs_full_res, downsample_local);
    }
    else if (downsample_local == 0)
      seqs = seqs_full_res;
    
    // Make sure that the images are resized to the correct resolution.
    disps = resize_images(disps, seqs[0].image.w, seqs[0].image.h);
  }

  // Convert the depth to disparity if we are working in depth-mode.
  if (!UseDisparity) {
    for (int i = 0; i < (int)disps.size(); ++i) {
      disps[i] = depth_to_disparity(disps[i]);
    }
  }

  
  // Write the output (full resolution).
  if(odisp.size() || odepth.size())
    write_disparities_and_depths(disps, odisp, odepth, seqs_full_res[0].image.w, seqs_full_res[0].image.h);
  
  // Reset the downsampled images.
  seqs = SeqType::downsample(seqs_full_res, downsample);
  disps = resize_images(disps, seqs[0].image.w, seqs[0].image.h);



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
  double flowImageCostWeight = kSceneFlowImageCostWeight;
  double dispConstWeight = 0.0;
  double dispImageCostWeight = kDispImageCostWeight;
  double flowBetaWeight = VariationalSceneFlowProblem<true>::defaultBeta.y;
  // double dispImageCostWeight = ???;
  
  std::string depthName = "", dispName = "", flowName = "";
  std::string debugRoot = "";
  std::string odisp = "", odepth = "";
  std::string oflow = "";

  int nclosest = -1;
  int imageIndex = 0;
  int useDepth = 0;

  cline.registerOption("downsample", "Number of times to perform downsampling.", &downsample, 'd');
  cline.registerOption("vdisp_dsmax", "Maximum downsample to reach for disparity.", &vdisp_max, 0);
  cline.registerOption("vdisp", "Do variational disparity", &variationalDisparity, 0);
  cline.registerOption("vflow", "Do variational flow", &variationalFlow, 0);
  cline.registerOption("use_truth", "Use truth", &useTruth, 0);

  cline.registerOption("vdisp_const", "Variational disparity constraint.", &dispConstWeight);
  cline.registerOption("dimg", "Image weight for disparity estimation.", &dispImageCostWeight);

  cline.registerOption("flow", "Flow input name", &flowName, 0);
  cline.registerOption("disp", "Name for isparity maps.", &dispName, 0);
  cline.registerOption("depth", "Name for input depth maps (instead of disparity).", &depthName, 0);

  cline.registerOption("index", "Image index to load", &imageIndex, 0);
  cline.registerOption("root", "Debug root", &debugRoot, 0);
  
  cline.registerOption("odisp", "Output name for disparity maps.", &odisp, 0);
  cline.registerOption("odepth", "Output name for depth maps.", &odepth, 0);
  cline.registerOption("oflow", "Output name for the flow maps.", &oflow, 0);
  
  cline.registerOption("fimg", "Image weight", &flowImageCostWeight, 0);
  cline.registerOption("fbeta", "Flow beta weight", &flowBetaWeight, 0);
  
  cline.registerOption("nclosest", "Number of images to use surrounding.", &nclosest, 0);
  cline.registerOption("use_depth", "Use depth formulation instead of disparity.", &useDepth, 0);

  cline.parse(ac, av);

  std::vector<SeqType> seqs_full_res = DisparitySequence::loadSequences(av + optind, ac - optind);
  std::vector<SeqType> seqs;

  if(imageIndex){
    for(int i=0; i<(int)seqs_full_res.size(); i++)
      seqs_full_res[i].load(imageIndex);
  }
    
  
  if(downsample)
    seqs = DisparitySequence::downsample(seqs_full_res, downsample);
  else
    seqs = seqs_full_res;

  if(nclosest < 0)
    nclosest = seqs.size();
  
  nclosest = std::max(1, nclosest);

  std::vector<Imagef> disps;
  
  if(dispName.size()){    
    printf("Using disps\n");

    for(int i=0; i<(int)seqs.size(); i++) {
      disps.push_back(nacb::Imagef((boost::format(dispName) % i).str().c_str()));
      disps.back() = disps.back().resize(seqs[i].image.w, seqs[i].image.h);
    }
  }
  else if(depthName.size()){
    // Load in the depth maps 
    printf("Using depth maps.\n");

    for(int i=0; i<(int)seqs.size(); i++) {
      Imagef depth = nacb::Imagef((boost::format(depthName) % i).str().c_str());
      std::cout << "Depth:" << depth.w << " " << depth.h << "\n";
      disps.push_back(disparity_to_depth(depth));
      disps.back() = disps.back().resize(seqs[i].image.w, seqs[i].image.h);
    }
  }
  else {
    fprintf(stderr, "Need starting disparity --disp\n");
    return 0;
  }

  // Load in the ground truth disparities.
  if(useTruth){
    printf("Loading ground truth disps.\n");
    disps.clear();
    
    for(int i=0; i<(int)seqs.size(); i++)
      disps.push_back(seqs[i].disp);
  }

  // Replace invisible regions with the mean.
  for (int i = 0; i < (int)disps.size(); ++i) {
    double mean = 0;
    int count = 0;
    for (int y = 0; y < seqs[i].image.h; ++y) {
      for (int x = 0; x < seqs[i].image.w; ++x) {
	if (seqs[i].image(x, y, 3) > 0) {
	  mean += disps[i](x, y);
	  count++;
	}
      }
    }
    mean /= count;
    
    for (int y = 0; y < disps[i].h; ++y) {
      for (int x = 0; x < disps[i].w; ++x) {
	if (seqs[i].image(x, y, 3) != 255){
	  float w = float(seqs[i].image(x, y, 3))/255.0f;
	  disps[i](x, y) = mean*(1.0 - w) + disps[i](x, y)*w;
	}
      }
    }

    printf("mean %d is: %f\n", i, mean);
  }

  // Write out the ground truth (as depth).
  for(int i=0; i<(int)seqs.size(); i++){
    const std::string fname = (boost::format("/tmp/input-%04d.rfi") % i).str();
    const int upw = seqs[i].disp.w << downsample;
    const int uph = seqs[i].disp.h << downsample;
    
    disparity_to_depth(seqs[i].disp).resize(upw, uph).write(fname.c_str());
  }

  if(debugRoot.size())
    bfs::create_directory(debugRoot);
  else 
    std::cout << "Not writing debug information.\n";

  // All things variational 
  {
    if(variationalDisparity){
      if (useDepth) {
	variational_disparity_main<false>(seqs_full_res, seqs, disps,
					 dispImageCostWeight, dispConstWeight,
					 variationalDisparity, nclosest,
					 downsample, vdisp_max, 
					 odisp, odepth, debugRoot);
      }
      else {
	variational_disparity_main<true>(seqs_full_res, seqs, disps,
					 dispImageCostWeight, dispConstWeight,
					 variationalDisparity, nclosest,
					 downsample, vdisp_max, 
					 odisp, odepth, debugRoot);
      }
    }
    
    if(variationalFlow){
      std::vector<Imagef> flows;
            
      for(int i=0; i<(int)seqs.size(); i++){
	nacb::Imagef dsource = disps[i];

	if (useDepth) {
	  std::cout << "Using depth.  Converting input disps to depth.  All flows will be depth representation.\n";
	  dsource = depth_to_disparity(dsource);
	}

	Imagef sf(dsource.w, dsource.h, 4);
	sf = 0;
	sf.setChannel(0, dsource);

	if(flowName.size()){
	  Imagef flowUse((boost::format(flowName) % i).str().c_str());
	  
	  flowUse = flowUse.resize(sf.w, sf.h);
	  
	  for(int k=0; k<3; k++)
	    sf.setChannel(k + 1, flowUse.getChannel(k));
	}
	flows.push_back(sf);
      }

      while(downsample >= vdisp_max){ 
	// The disparity-based implementation.
	if (!useDepth) {
	  scene_flow_output<true>::func_t output_func = 
	    boost::bind(::variational_scene_flow_output<DisparityOffsetBase<true> >, debugRoot, _1, seqs, _2, _3, _4, _5, _6);
	  
	  for(int i=0; i<variationalFlow; i++)
	    variational_scene_flow<DisparityOffsetBase<true> >(seqs, flows, downsample, 
							       flowImageCostWeight*sqrt(double(1 << downsample)), 
							       flowBetaWeight, dispConstWeight, output_func, nclosest);
	}
	else {
	  // The depth-based implementation.
	  scene_flow_output<false>::func_t output_func = 
	    boost::bind(::variational_scene_flow_output<DisparityOffsetBase<false> >, debugRoot, _1, seqs, _2, _3, _4, _5, _6);
	  
	  for(int i=0; i<variationalFlow; i++)
	    variational_scene_flow<DisparityOffsetBase<false> >(seqs, flows, downsample, 
								flowImageCostWeight*sqrt(double(1 << downsample)), 
								flowBetaWeight, dispConstWeight, output_func, nclosest);

	}
	downsample--;
	std::cout << "New downsample: " << downsample << "\n";

	if(downsample > 0){
	  seqs = SeqType::downsample(seqs_full_res, downsample);
	}
	else if (downsample == 0)
	  seqs = seqs_full_res;
	
	// Make sure that the images are resized to the correct resolution.
	flows = resize_images(flows, seqs[0].image.w, seqs[0].image.h);
      }

      if (downsample != 0) {
	seqs = seqs_full_res;
	flows = resize_images(flows, seqs[0].image.w, seqs[0].image.h);
      }

      save_depth_and_flow(flows, odepth, oflow, useDepth == 0);
    }
  }
  return 0;
}
