#include <GL/glew.h>
#include "utils.h"
#include "VariationalBasisProblemUV.h"
#include "DisplaceUV.h"
#include "DisparitySequence.h"
#include "DisplaceOffsetBasis.h"

#include "utigl/glwindow.h"
#include "utigl/fbo.h"
#include "ik/mesh.h"
#include "ik/triangle.h"
#include <nimage/image.h>

#include <nmisc/commandline.h>

using namespace nacb;
using namespace std;

/**
   A set of options for the multi-res basis flow.
*/
class MultiResOptions {
public:
  MultiResOptions() {
    texw = 512;
    texh = 512;

    imageCostFlowWeight = 1.0;
    imageCostStereoWeight = 1.0;
    flowReg = 0.5;
    stereoReg = 0.5;

    downsample = 0;
    downsampleMax = 0;
    maxits = 4;

    stereoTerms = 1;
    referenceTime = 0;
    index1 = 0;
    index2 = 7;

    reweightPower = 0.5;
  }

  double reweightPower;
  double imageCostFlowWeight;
  double imageCostStereoWeight;
  double flowReg;
  double stereoReg;
  
  int downsample;
  int downsampleMax;
  int maxits;

  int referenceTime;
  int index1;
  int index2;
  
  int texw, texh;
  int stereoTerms;

};


using namespace DisplaceOffsetBasis;


nacb::Imagef multiResolutionBasisFlow(std::vector<DisparitySequence>& seqs_full_res,
				      FlowBasis3D::ptr& basis,
				      Mesh& mesh, 
				      BaseMeshAnimation& meshAnimation,
				      MultiResOptions& opts,
				      nacb::Imagef& depthOffset) {
  typedef ImageBaseSkinnedGeometry BaseGeometryType;
  DisplacedMesh dispMesh;
  
  typedef VariationalBasisProblemUV VariationalBasisProblemType;
  std::vector<DisparitySequence> seqs;

  depthOffset = depthOffset.resize(depthOffset.w >> opts.downsample,
				   depthOffset.h >> opts.downsample);				     
  
  if(opts.downsample)
    seqs = DisparitySequence::downsample(seqs_full_res, opts.downsample);
  else
    seqs = seqs_full_res;

  std::vector<std::vector<nacb::Imagef> >  images;
  std::vector<DisplaceOffsetUV>     warps;
  std::vector<std::vector<nacb::Imagef> >    weights;
  std::vector<pair<int, int> >   pairs;
  


  // Initialize all stereo pairs.
  if (opts.stereoTerms) {
    pairs = buildStereoPairs(seqs.size());
  }

  std::vector<int> times;
  for (int t=opts.index1; t <= opts.index2; t++) {
    times.push_back(t);
  }

  int downsample_local = opts.downsample;
  while (downsample_local >= opts.downsampleMax) {    
    BaseGeometry::ptr baseGeometrySubdiv = 
      createBaseGeometryFromMesh<BaseGeometryType>(mesh, 
				 opts.texw >> downsample_local, 
				 opts.texh >> downsample_local, dispMesh);
    depthOffset = depthOffset.resize(baseGeometrySubdiv->getWidth(),
				     baseGeometrySubdiv->getHeight());

    meshAnimation.setBaseGeometry(baseGeometrySubdiv);

    initImagesWarpsAndWeights(seqs, basis, meshAnimation, dispMesh, baseGeometrySubdiv, times, depthOffset,
			      images, warps, weights, true);
    
    printf("Solving problem.\n");
    // Solve the problem
    Vec2f imageCostWeightPairDS(opts.imageCostStereoWeight * pow(double(1 << downsample_local), opts.reweightPower) / pairs.size(),
				opts.imageCostFlowWeight * pow(double(1 << downsample_local), opts.reweightPower) / images.size());
    Vec2f regWeight(opts.stereoReg, opts.flowReg);

    for (int k=0; k<opts.maxits; k++) {
      //depthOffset = depthOffset.medianFilter(4);

      VariationalBasisProblemType vari(warps, images, basis, 
				       opts.referenceTime, meshAnimation , weights, 
				       pairs, depthOffset, times, imageCostWeightPairDS, regWeight);
      nacb::Imagef residuals;
      mg::wcycle_solver solver(2, 1, 4, 4, true);
      
      
      printf("%dx%d  %lf %lf %lf\n", depthOffset.w, depthOffset.h,
	     VariationalBasisProblemType::imageDataEnergy(depthOffset, images, warps, weights, pairs, 
							  times, opts.referenceTime, basis, meshAnimation),
	     VariationalBasisProblemType::smoothnessEnergy(depthOffset),
	     VariationalBasisProblemType::smoothnessEnergyFlow(depthOffset));
	    
	     
      
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
	std::vector<nacb::Imagef> warpedImages;

	for (int t=0; t<(int)times.size(); t++) {
	  std::string name = (boost::format("/tmp/res-%04d-%d.png") % k % t).str();
	  warpedImages.push_back(warpBack(images[k], basis, meshAnimation, warps[k], depthOffset, t));
	  warpedImages.back().write(name.c_str());
	}

	/*
	image_stats_t stats;
	for (int j=0; j<int(times.size()-1); j++) {
	  stats += image_compare(weights[k][j] * weights[k][j+1], warpedImages[j], warpedImages[j+1]);
	}
	std::cout << "seq" << k << "(flow):" << stats << std::endl;
	*/
      }
      
      printf("%dx%d %lf %lf %lf\n", depthOffset.w, depthOffset.h,
	     VariationalBasisProblemType::imageDataEnergy(depthOffset, images, warps, weights, pairs,
							  times, opts.referenceTime, basis, meshAnimation),
	     VariationalBasisProblemType::smoothnessEnergy(depthOffset),
	     VariationalBasisProblemType::smoothnessEnergyFlow(depthOffset));

      basis->save("/tmp/basis_results.b3d", depthOffset);
    }
    downsample_local--;
    if (downsample_local > 0)
      seqs = DisparitySequence::downsample(seqs_full_res, downsample_local);
    else if (downsample_local == 0)
      seqs = seqs_full_res;
  }

  BaseGeometry::ptr baseGeometrySubdiv = 
    createBaseGeometryFromMesh<BaseGeometryType>(mesh, 
						 opts.texw >> downsample_local, 
						 opts.texh >> downsample_local, dispMesh);

  meshAnimation.setBaseGeometry(baseGeometrySubdiv);

  depthOffset = depthOffset.resize(baseGeometrySubdiv->getWidth(),
				   baseGeometrySubdiv->getHeight());

  for (int t=0; t<(int)times.size(); t++) {
    nacb::Imagef dispMap = basis->getDisplacementMap(depthOffset, t);
    dispMesh.displaceAndOffset(baseGeometrySubdiv, dispMap);
    dispMesh.saveobj((boost::format("/tmp/displaced-%d.obj") % t).str().c_str(), false, true);    
  }

  return depthOffset;
}


/** \brief Quick function for testing whether the basis function looks correct.
 */
void showBasis(const FlowBasis3D::ptr & basis, int numTimes) {
  // Evaluate the basis functions at several times.
  nacb::Imagef D(1, 1, basis->getNumDimensions());
  nacb::Matrix basisFunctions[4] = {nacb::Matrix::zeros(basis->getNumDimensions(), numTimes),
				    nacb::Matrix::zeros(basis->getNumDimensions(), numTimes),
				    nacb::Matrix::zeros(basis->getNumDimensions(), numTimes),
				    nacb::Matrix::zeros(basis->getNumDimensions(), numTimes)};
  D = 0;

  for (int i=0; i<basisFunctions[0].m; i++) {
    for (int t=0; t<basisFunctions[0].n; t++) {
      D(0, 0, i) = 1;
      double displace = 0;
      nacb::Vec3f offset = (*basis)(displace, 0, 0, D, t);
      D(0, 0, i) = 0;
      
      basisFunctions[0](i, t) = displace;
      basisFunctions[1](i, t) = offset.x;
      basisFunctions[2](i, t) = offset.y;
      basisFunctions[3](i, t) = offset.z;
    }
  }
  basisFunctions[0].printMatlab("M0");
  basisFunctions[1].printMatlab("M1");
  basisFunctions[2].printMatlab("M2");
  basisFunctions[3].printMatlab("M3");
}


int main(int ac, char * av[]){
  std::string geomFile;
  std::string basisFile;
  std::string outputFile;

  nacb::CommandLine cline;
  MultiResOptions opts;
  std::string basisType = "standard";

  std::string bonesFileName;
  std::string animFileName;
  int numBasis = 0;
  int useDepthBasis = 0;

  cline.registerOption("geom", "The geometry file", &geomFile, 'g');
  cline.registerOption("img", "Disparity image weight (stereo term, will be normalized by pairs)", &(opts.imageCostStereoWeight), 'i');
  cline.registerOption("flow", "Flow image weight (will be normalized by number of seqs)", &(opts.imageCostFlowWeight), 'f');

  cline.registerOption("downsample", "Downsample value.", &(opts.downsample), 'd');
  cline.registerOption("vdisp_dsmax", "Downsample max.", &(opts.downsampleMax), 0);
  cline.registerOption("reference", "Reference time index.", &(opts.referenceTime), 0);
  cline.registerOption("index1", "Image index 1.", &(opts.index1), 0);
  cline.registerOption("index2", "Image index 2.", &(opts.index2), 0);
  cline.registerOption("maxits", "Maximum iterations.", &(opts.maxits), 0);
  cline.registerOption("basisFile", "Initialize to this displacement.", &basisFile, 0);
  cline.registerOption("output", "Output displacement name.", &outputFile, 0);
  cline.registerOption("flowReg", "Flow regularization", &(opts.flowReg), 0);
  cline.registerOption("stereoReg", "Stereo regularization", &(opts.stereoReg), 0);
  cline.registerOption("numBasis", "Number of basis elements", &numBasis, 0);
  cline.registerOption("basis", "Basis type (cosine, standard)", &basisType, 0);
  cline.registerOption("bones", "Bones file name", &bonesFileName, 0);
  cline.registerOption("anim", "Animation of bones", &animFileName, 0);
  cline.registerOption("stereo", "Use the stereo terms", &(opts.stereoTerms), 0);
  cline.registerOption("useDepthBasis", "Use depth basis (cosine offset only)", &useDepthBasis, 0);
  cline.registerOption("texw", "Texture width", &(opts.texw));
  cline.registerOption("texh", "Texture height", &(opts.texh));
  cline.registerOption("reweight", "Reweight power.", &(opts.reweightPower));
  cline.parse(ac, av);

  if (optind >= ac) {
    std::cout << "Need at least two image sequences.\n";
    return 1;
  }

  if (!geomFile.size()) {
    std::cout << "Need a geometry file.\n";
    return 2;
  }

  Mesh mesh;
  mesh.loadobj(geomFile.c_str());

  // Load in the animation (and bones)
  BaseMeshAnimation baseMeshAnimation;
  if (animFileName.size() && bonesFileName.size()) {
    PoseKeys poseKeys;
    baseMeshAnimation = BaseMeshAnimation(mesh);

    if (!poseKeys.load(animFileName.c_str())) {
      std::cerr << "Error loading animation from:" << animFileName << std::endl;
      return 0;
    }
    baseMeshAnimation.setPoseKeys(poseKeys);

    if (!baseMeshAnimation.loadArmature(bonesFileName.c_str())) {
      std::cerr << "Error loading armature from " << bonesFileName << std::endl;
      return 0;      
    }
  }

  printf("%ld %ld %ld\n", mesh.vert.size(), mesh.tris.size(), mesh.tvert.size());

  std::vector<DisparitySequence> seqs_full_res = 
    DisparitySequence::loadSequences(av + optind, ac - optind);
 

  GLWindow * window = new GLWindow(64, 64);//opts.texw, opts.texh);
  glewInit();
  FrameBufferObject fbo(max(1024, opts.texw), max(opts.texh, 1024));
  fbo.bind(1);

  nacb::Imagef result;
  FlowBasis3D::ptr basis(new StandardFlowBasis3D(opts.referenceTime));

  if (basisType == "cosine") {
    basis = FlowBasis3D::ptr(new CosineFlowBasis3D(opts.referenceTime, numBasis, (opts.index2 - opts.index1 + 1), useDepthBasis));
  }
  else if (basisType == "perframe") {
    basis = FlowBasis3D::ptr(new PerFrameFlowBasis3D(opts.referenceTime, (opts.index2 - opts.index1 + 1)));
  }
  else if (basisType == "cosinedepth") {
    basis = FlowBasis3D::ptr(new CosineDepthOffset(opts.referenceTime, numBasis, (opts.index2 - opts.index1 + 1)));
  }
  else if (basisType == "standard2D") {
    basis = FlowBasis3D::ptr(new StandardFlowBasis2D(opts.referenceTime));
  }
      
  if (basisFile != "") {
    nacb::Imagef resultInit;
    FlowBasis3D::ptr basisInit = FlowBasis3D::create(basisFile, resultInit);

    if (!basisInit) {
      std::cout << "Cannot load result from:" <<  basisFile << std::endl;
      return -1;
    }

    showBasis(basisInit, (opts.index2 - opts.index1 + 1));

    // I used to just enforce the same initialization as before, but now
    // it projects the input the desired basis (much better for bootstrapping).
    //
    // Use the first basis as initialization.
    std::vector<nacb::Imagef> maps;
    std::vector<int> times;
    
    // Extract the displacement maps.
    for (int t = opts.index1; t <= opts.index2; t++) {
      maps.push_back(basisInit->getDisplacementMap(resultInit, t));
      times.push_back(t);
    }
    
    nacb::Matrix B = basis->getBasisMatrix(times).transpose();
    nacb::Matrix coeff(4*times.size(), 1);

    result = nacb::Imagef(resultInit.w, resultInit.h, basis->getNumDimensions());

    // Project onto the basis.
    for (int y=0; y<maps[0].h; y++) {
      for (int x=0; x<maps[0].w; x++) {
	for (int t=0; t<(int)times.size(); t++) {
	  for (int k=0; k<4; k++) 
	    coeff[4*t + k] = maps[t](x, y, k);
	}
	//printf("%d x %d  %d x %d  %d x %d\n", x, y, result.w, result.h, maps[0].w, maps[0].h);
	nacb::Matrix M = nacb::Matrix::LlinLeastSq(B, coeff);
	//M.printMatlab("M");
	//printf("%d %d\n", basis->getNumDimensions(), result.nchannels);
	
	for (int k=0; k<M.m; k++)
	  result(x, y, k) = M[k];
      }
    }
    

  }
  else {
    result = nacb::Imagef(opts.texw, opts.texh, basis->getNumDimensions());
  }

  result = multiResolutionBasisFlow(seqs_full_res, basis, mesh, baseMeshAnimation, opts, result);

  if (outputFile.size()) {
    result.resize(opts.texw, opts.texh).write(outputFile.c_str());
    result.resize(opts.texw, opts.texh).getChannel(0).getNormalizedImage().write("/tmp/d-final.png");
    result.resize(opts.texw, opts.texh).getChannel(2).getNormalizedImage().write("/tmp/dv-final.png");
  }

  delete window;
  return 0;
}
