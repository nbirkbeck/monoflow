/**
   This is a small program to test the approximation of 
   a dpl displacement file with a cosine basis.

   (c) Neil Birkbeck, Feb 2011
*/

#include "ik/displace.h"
#include <stdexcept>
#include "ik/mesh.h"
#include "ik/defmodel.h"
#include "ik/armature.h"
#include "ik/MeshMotionConstraint.h"
#include <nmisc/commandline.h>
#include <nmath/vec3.h>
#include <nmath/matrix.h>
#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"
#include "autorecon/clbfile.h"
#include "FlowBasis3D.h"


using namespace dpl;


/** \brief Return factorial of x (x!)
 */
double fact(int x) {
  if (x <= 1) return 1.0;

  double val = 1.0;
  for (int i=1; i<=x; i++) {
    val *= i;
  }
  return val;
}


/** \brief return x choose y
 */
double binomial(int x, int y) {			
  if (x <= 1) return 1.0;
  if (y < x/2) return binomial(x, x - y);

  double num = 1.0;
  for (double v=y; v<=x; v++) {
    num *= v;
  }
  return num/fact(x - y);
}


/** \brief Compute the bezier basis (directly).
 */
nacb::Matrix getBezierBasisMatrixSimple(int numTimes, int numBasis) {
  nacb::Matrix B = nacb::Matrix::zeros(numTimes, numBasis);

  // use the appropriate n (e.g., the first n, s.t.
  int n = numBasis - 1;
    
  for (int i=0; i<numBasis; i++) {
    int k = i;
    double weight = binomial(n, k);
    // printf("weight: %f  %d, %d  %f\n", weight, n, k, fact(n)/(fact(k)*fact(n - k)));
    
    for (int j=0; j<numTimes; j++) {
      double t = double(j)/double(numTimes - 1);
      B(j, i) = weight * pow(t, k) * pow((1.0 - t), n - k);
    }   
  }
  return B;
}



/** Compute the bezier basis matrix.  This is ill-conditioned for a large number
    of basis elements.

    This is the dejastro algorithm (or whatever his name is).  Trying to 
    get a better approximation, but I think the basis matrix is just ill-conditioned.
*/
nacb::Matrix getBezierBasisMatrix(int numTimes, int numBasis) {
  nacb::Matrix B = nacb::Matrix::zeros(numTimes, numBasis);

  // use the appropriate n (e.g., the first n, s.t.
  int n = numBasis - 1;
    
  for (int i=0; i<numBasis; i++) {
    std::vector<double> coeffs(n + 1, 0); 
    coeffs[i] = 1;

    for (int j=0; j<numTimes; j++) {
      double t0 = double(j)/double(numTimes - 1);

      std::vector<double> coeffsUse = coeffs;

      // Eval at t0 (numerically stable)
      while (coeffsUse.size() != 1) {
	std::vector<double> newCoeffs(coeffsUse.size() - 1, 0);
	for (int ii=0; ii<newCoeffs.size(); ii++) {
	  newCoeffs[ii] = (1.0 - t0) * coeffsUse[ii] + t0 * coeffsUse[ii + 1];
	}
	coeffsUse = newCoeffs;
      }

      B(j, i) = coeffsUse[0];
    }   
  }
  return B;
}


/** \brief Compute the cosine basis matrix.
 */
nacb::Matrix getCosineBasisMatrix(int numTimes, int numBasis) {
  CosineBasis basis(numTimes, 0);
  nacb::Matrix B(numTimes, numBasis);
  for (int i=0; i<numBasis; i++)
    for (int j=0; j<numTimes; j++) 
      B(j, i) = basis.evaluateBasis(i, j);
  return B;
}


/** \brief Compute the catmull rom basis matrix.
    
    The basis coefficients are really the control points.  These
    are spaced such that the required number of basis elements 
    fits up to the total number of times.
 */
nacb::Matrix getCatmullRomBasisMatrix(int numTimes, int numBasis) {
  double frameSpacing = (double(numTimes - 1)/ (numBasis - 1));
  printf("spacing: %f  max: %d (%d)\n", frameSpacing, int(frameSpacing * (numBasis - 1)), numTimes);
  nacb::Matrix B = nacb::Matrix::zeros(numTimes, numBasis);
  
  // Catmull rom is only dependent on the control points.  Uses tangents
  // and control points.  From this we can compute the basis elements of a point,
  // as it affects either the point directly or the tangent. This is done for 
  // intervals.
  for (int i=0; i<numBasis; i++) {
    for (int j = max(0, (int)floor(frameSpacing*i - 2*frameSpacing));
	 j<min(numTimes, (int)ceil(frameSpacing*i + 2*frameSpacing)); j++) 
    {
      double intervald = (double(j)/frameSpacing);
      int interval = (int)intervald;
      double t = intervald - interval;
      // p0, m0, p1, m1
      double coeffs[4] = {(2.0*t*t*t - 3*t*t + 1),
			  t*t*t - 2.0*t*t + t,
			  -2.0*t*t*t + 3*t*t,
			  t*t*t - t*t};
      // m0 = p_{t+1} - p_{t-1}
      // m1 = p_{t+2} - p_{t}
      if (interval == i - 2) {
	// affects m1
	B(j, i) = coeffs[3]/frameSpacing;
      }
      else if (interval == i - 1) {
	// Affects p1 and m0 
	B(j, i) = coeffs[2] + coeffs[1]/frameSpacing;
      }
      else if (interval == i) {
	// Affects p0 and m1
	B(j, i) = coeffs[0] - coeffs[3]/frameSpacing;
      }
      else if (interval == i + 1) {
	// affects m0 
	B(j, i) = - coeffs[1]/frameSpacing;
      }
    }
  }
  
  return B;
}


/** \brief Approximate the input displacements with the given basis.
    
    \return  The approximated displacements (basis coefficients are thrown out).
 */
std::vector<dpl::displace_t> approximateDisplacements(const std::vector<dpl::displace_t> & disp,
						      nacb::Matrix& B) {
  std::vector<dpl::displace_t> results = disp;
  //nacb::Matrix B(disp.size(), B.n);
  nacb::Matrix R = B*(B.transpose()*B).inverse()*B.transpose();

  for (int pointIndex=0; pointIndex < disp[0].size(); pointIndex++) {
    nacb::Matrix X[3] = {nacb::Matrix(disp.size(), 1),
			 nacb::Matrix(disp.size(), 1),
			 nacb::Matrix(disp.size(), 1)};


    //if (pointIndex % 1000 == 0) std::cout << pointIndex << std::endl;
    // Grab all the observations (x, y, z)
    for (int t=0; t<disp.size(); t++) {
      for (int k=0; k<3; k++) {
	X[k][t] = disp[t][pointIndex].data[k]; 
      }
    }

    nacb::Matrix Y[3] = {nacb::Matrix(1, 1),
			 nacb::Matrix(1, 1),
			 nacb::Matrix(1, 1)};
    
    // Reconstruc each component (x, y, z) (project onto basis, and back out)
    // B = X
    // B'*B = B'*X
    // C = inv(B'*B)*(B*X')
    // Y = B*C = (B*inv(B'*B))*(B*X)
    for (int k=0; k<3; k++) {
      //Y[k] = B * nacb::Matrix::LlinLeastSq(B, X[k]);
      Y[k] = R * X[k]; // Much faster than solving at each frame.
    }
    
    // Unpack the results into the displacements.
    for (int t=0; t<disp.size(); t++) {
      for (int k=0; k<3; k++) {
	results[t][pointIndex].data[k] = Y[k][t];
      }
    }    
  }
  return results;
}


/** \brief Compute the SSD distance between two displacement sequences.
 */
double dplDistance(const std::vector<displace_t> & d1,
		   const std::vector<displace_t> & d2)
{
  if (d1.size() != d2.size()) {
    throw std::runtime_error("Displacements of different lengths cannot be compared.");
  }

  double error = 0.0;
  int cnt = 0;

  for (int i=0; i<d1.size(); i++) {
    if (d1[i].size() != d2[i].size()) {
      throw std::runtime_error("Displacements of different lengths cannot be compared.");
    }

    for (int j=0; j<d1[i].size(); j++) {
      error += (d1[i][j] - d2[i][j]).len();
      cnt++;
    }
  }
  return error/cnt;
}


/** \brief This function tries to compensate for any motion in the displacements
           by using the kinematic motion.
   
   Update the animation and the ground truth displacements so that the "recon"
   displaecements best fit the current displaced object.

   This does in fact reduce the residuals mostly when using fewer basis functions.
*/
PoseKeys updateAnimationAndDisplacement(Mesh& mesh, 
					Armature& arm,
					PoseKeys& keys,
					std::vector<dpl::displace_t>& truth,
					std::vector<dpl::displace_t>& recon)
{
  PoseKeys result(arm.getBoneNames());
  printf("Num frames: %d\n", truth.size());

  for (int frame=0; frame<truth.size(); frame++) {
    printf("%d ", frame);
    fflush(stdout);

    std::vector<nacb::Vec3f> restVert = mesh.restVert;
    std::vector<nacb::Vec3f> vert = mesh.vert;

    // Get the current ground truth displaced mesh.
    {
      // Update the rest vert to the current ground truth.
      for (int k=0; k<mesh.restVert.size(); k++) {
	mesh.restVert[k] = restVert[k] + truth[frame][k];
      }
      
      // Animate the mesh;
      keys.setCurrent(frame);
      arm.set(keys);
      arm.animate(mesh);
    }

    // Save the displaced vertices.
    std::vector<nacb::Vec3f> displacedVert = mesh.vert;
    mesh.vert = vert;

    // Fit the animation using the reconstructed rest vert.
    for (int k=0; k<mesh.restVert.size(); k++) {
      mesh.restVert[k] = restVert[k] + recon[frame][k];
    }

    MeshMotionConstraintSolver solver(arm, mesh, displacedVert);

    // Restore the rest vert (needed to make sure that displacements get extracted).
    mesh.restVert = restVert;
    mesh.vert = displacedVert; 

    // Get the new ground truth displacements using the updated animation.
    truth[frame] = extractFromPosed(arm, mesh);

    // Restore the mesh vertices.
    mesh.vert = vert;

    result.append(&arm);
  }
  printf("\n");
  return result;
}


Vec3f mult(const nacb::Matrix& P, const nacb::Vec3f & X) {
  nacb::Vec3f PX(0, 0, 0);
  for (int i=0; i<3; i++) {
    PX.data[i] = P(i, 3);

    for (int j=0; j<3; j++) {
      PX.data[i] += P(i, j)*X.data[j];
    }
  }
  return PX;
}


std::vector<nacb::Vec3f>  linearSolveDisplacement(nacb::Matrix & B,
						  std::vector<nacb::Matrix>& Ps,
						  const std::vector<nacb::Vec2f> & proj,
						  std::vector<nacb::Vec3f> & Xs,
						  const std::vector<int> & times) {
  // | \Pi(P * T * (X + d)) - [x;y] |^2  
  // 1) px - x * pz = 0
  // 2) py - y * pz = 0
  // 
  // px = (P * T)_1 * (X + d)
  // 
  // 1) (P * T)_1 * d - x * (P * T)_3 * d = x*(P * T)_3 * X - (P * T)_1 * X
  // 2) (P * T)_2 * d - x * (P * T)_3 * d = x*(P * T)_3 * X - (P * T)_2 * X
  // d = \sum lamba_i B_i(t)
  
  nacb::Matrix A = nacb::Matrix::zeros(2*Ps.size(), B.n*3);
  nacb::Matrix b = nacb::Matrix::zeros(2*Ps.size(), 1);
  
  // Build constraints for each time step.
  for (int t=0; t<Ps.size(); t++) {
    nacb::Vec3f PX = mult(Ps[t], Xs[t]);
    double x = proj[t].x, y = proj[t].y;
    //printf("%f %f  %f %f\n", x, y, PX.x/PX.z, PX.y/PX.z);
    for (int k=0; k<B.n; k++) {
      A(2*t, 3*k + 0) = Ps[t](0, 0)*B(times[t], k) - x*Ps[t](2, 0)*B(times[t], k) ;
      A(2*t, 3*k + 1) = Ps[t](0, 1)*B(times[t], k) - x*Ps[t](2, 1)*B(times[t], k) ;
      A(2*t, 3*k + 2) = Ps[t](0, 2)*B(times[t], k) - x*Ps[t](2, 2)*B(times[t], k) ;

      b[2*t] = x*PX.z - PX.x; // The following component is already in rhs: + (x*Ps[t](2, 3) - Ps[t](0, 3));

      
      A(2*t + 1, 3*k + 0) = Ps[t](1, 0)*B(times[t], k) - y*Ps[t](2, 0)*B(times[t], k) ;
      A(2*t + 1, 3*k + 1) = Ps[t](1, 1)*B(times[t], k) - y*Ps[t](2, 1)*B(times[t], k) ;
      A(2*t + 1, 3*k + 2) = Ps[t](1, 2)*B(times[t], k) - y*Ps[t](2, 2)*B(times[t], k) ;

      b[2*t + 1] = y*PX.z - PX.y; // The following component is already in rhs: + (y*Ps[t](2, 3) - Ps[t](1, 3));
    }
  }
  //A.printMatlab("A");
  // Get the basis coefficients.
  nacb::Matrix C = nacb::Matrix::LlinLeastSq(A, b);
  
  nacb::Matrix R = A*C - b;
  //printf("res: %f (was: %f)\n", sqrt(R.dot(R)/Ps.size()), sqrt(b.dot(b)/Ps.size()));

  // Reconstruct the points from the basis coefficients.
  std::vector<nacb::Vec3f> disp(Xs.size());

  // Reshape the coefficients to be (m x 3)
  C.n = 3;
  C.m = B.n;


  nacb::Matrix X = B * C;

  for (int t=0; t<Ps.size(); t++) {
    //printf("%f %f %f\n", X(t, 0), X(t, 1), X(t, 2));
    disp[t] = nacb::Vec3f(X(t, 0), X(t, 1), X(t, 2));
  }

  return disp;
}



/** \brief Use the basis functions in the reconstruction.
 */
std::vector<dpl::displace_t>
reconstructDisplacements(const std::vector<nacb::Matrix> & P,
			 Armature &arm,
			 PoseKeys & keys,
			 Mesh& mesh, nacb::Matrix & B,
			 const std::vector<std::vector<std::vector<nacb::Vec2f> > > & projections,
			 double plast = 1.0) { 
  int numFrames = B.m;
  
  // For each point in the mesh, reconstruct the displacement
  // of the point from a sequence of observations in a single view.

  std::vector<std::vector<Mat4x4> > Ts;
  std::vector<int> times;

  for (int t=0; t<numFrames; t++) {
    keys.setCurrent(t);
    arm.set(keys);
    arm.animate(mesh);
    
    Ts.push_back(arm.getTransforms());
  }

  std::vector<dpl::displace_t> recon(numFrames, dpl::displace_t(mesh.vert.size()));
  int missing = 0;

  for (int i=0; i<mesh.vert.size(); i++) {
    std::vector<nacb::Matrix> Ps; //! The combined projection matrix.
    std::vector<Vec3f> Xs;        //! The transformed rest vert.
    std::vector<nacb::Vec2f> proj;

    for (int t=0; t<numFrames; t++) {
      // For each frame, compute the transform T.
      nacb::Mat4x4 T = nacb::Mat4x4::zeros();
      nacb::Vec3f X(0, 0, 0);
      double weight = 0;
      
      // Accumulate over all influencing bones.
      for (int k=0; k<mesh.bone_weights[i].size(); k++) {
	T = T + Ts[t][mesh.bone_weights[i][k].bone] * mesh.bone_weights[i][k].weight;
	weight += mesh.bone_weights[i][k].weight;
      }
      T = T*(1.0/weight);


      // The image observation is \PI(P * T * (X + d))
      std::vector<Matrix> xs;
      for (int cam = 0; cam < P.size(); cam++) {
	nacb::Matrix PT = P[cam] * T;

	if (cam == 0 || drand48() <= plast) 
	{
	  Xs.push_back(mesh.restVert[i]);
	  Ps.push_back(PT);
	  times.push_back(t);
	  
	  proj.push_back(projections[cam][i][t]);
	}
	else 	 
	  missing++;

      }

      /*
      // Test: get an approximation using triangulation per frame only.
      nacb::Vec3f v = triangulate(P[0], P[1],
				  projections[0][i][t],
				  projections[1][i][t]);
      recon[t][i] = dpl::itransform(Ts[t], mesh.bone_weights[i], v, 1) - mesh.restVert[i];
      */
    }
	 
    // Now get a linear solution to displacement.
    std::vector<nacb::Vec3f> r = linearSolveDisplacement(B, Ps, proj, Xs, times);
    for (int t=0; t<recon.size(); t++) {
      recon[t][i] = r[t];
    }
  }
  printf("Missing: %d\n", missing);
  return recon;
}


std::vector<std::vector<std::vector<nacb::Vec2f> > >
getImageProjections(const std::vector<nacb::Matrix> & P,
		    Armature &arm,
		    PoseKeys & keys,
		    Mesh& mesh,
		    const std::vector<dpl::displace_t>& disp,
		    double noiseAmount  = 0.0) {

  // For each time in the sequence, get image projection of all mesh points.
  int numFrames = std::min((int)disp.size(), (int)keys.getNumFrames());
  RawDeformationModel defmodel(disp);

  std::vector<std::vector<std::vector<nacb::Vec2f> > >allProjections;

  for (int cam=0; cam<P.size(); cam++) {
    std::vector<std::vector<nacb::Vec2f> > projections(mesh.vert.size(), 
						       std::vector<nacb::Vec2f>(numFrames));
    
    for (int t=0; t<numFrames; t++) {
      keys.setCurrent(t);
      arm.set(keys);
      
      defmodel.deform(mesh, arm, t);
      
      for (int i=0; i<mesh.vert.size(); i++) {
	projections[i][t] = project(P[cam], mesh.vert[i]);
	
	if (noiseAmount > 0) {
	  projections[i][t].x += (drand48() - 0.5)*noiseAmount;
	  projections[i][t].y += (drand48() - 0.5)*noiseAmount;
	}
      }
    }			 
    allProjections.push_back(projections);
  }
  return allProjections;
}


/** \brief Compute the indices that are attached to the bones.
 */
std::vector<int> getAttachedIndices(Mesh& mesh, std::set<int>& bones) {
  std::vector<int> indices;

  for (int i=0; i<mesh.vert.size(); i++) {
    double w = 0;
    for (int j=0; j<mesh.bone_weights[i].size(); j++) {
      if (bones.count(mesh.bone_weights[i][j].bone)) {
	w += mesh.bone_weights[i][j].weight;
      }
    }
    if (w > 0.97 && drand48() < 0.15) {
      indices.push_back(i);
    }
  }
  
  return indices;
}


/** \brief Return a vector with every "skip" indices.
 */
std::vector<int> every(Mesh& mesh, int skip) {
  std::vector<int> indices;

  for (int i=0; i<mesh.vert.size(); i+=skip) {
    indices.push_back(i);
  }
  
  return indices;
}


template <class T> T extractPart(const T&, std::vector<int>& indices);

template <>
std::vector<nacb::Vec3f> extractPart(const std::vector<nacb::Vec3f> & input, 
				     std::vector<int> &indices) {
  std::vector<nacb::Vec3f> result(indices.size());

  for (int i=0; i<indices.size(); i++) {
    result[i] = input[indices[i]];
  }
  return result;
}


template <>
std::vector<dpl::displace_t> extractPart(const std::vector<dpl::displace_t> & input, 
					 std::vector<int> &indices) {
  std::vector<dpl::displace_t> result = input;

  for (int i=0; i<input.size(); i++) {
    result[i] = extractPart(input[i], indices);
  }
  return result;
}


int main(int ac, char * av[])
{
  using namespace dpl;

  nacb::CommandLine cline;
  std::string geomFileName;
  std::string dplFile;
  std::string animFileName;
  std::string bonesFileName;
  std::string basisType = "cosine";
  std::string calibFileName;

  int maxNumBasis = -1;
  int minNumBasis = 1;
  int numIterations = 1;
  double noiseAmount = 0;
  double plast = 1.0;
  Mesh mesh;
  Armature arm;
  PoseKeys keys;

  cline.registerOption("plast", "The probability that the later cameras have observations", &plast, 0);
  cline.registerOption("numIterations", "Refine animation.", &numIterations, 0);
  cline.registerOption("anim", "The animation file.", &animFileName, 0);
  cline.registerOption("bones", "The bones file.", &bonesFileName, 0);
  cline.registerOption("geom", "The name of the geometry.", &geomFileName, 0);
  cline.registerOption("dpl", "The name of the displacement file.", &dplFile, 0);
  cline.registerOption("numBasis", "The maximum number of basis.", &maxNumBasis, 0);
  cline.registerOption("minBasis", "The min number of basis.", &minNumBasis, 0);
  cline.registerOption("basisType", "The basis type.", &basisType, 0);
  cline.registerOption("calib", "The image calibration file to use in synthetic image reconstruciton.", &calibFileName, 0);
  cline.registerOption("noise", "The amount of noise added to image observations.", &noiseAmount);
  cline.parse(ac, av);

  if (!dplFile.size()) {
    std::cerr << "Need displacements (--dpl)" << std::endl;
    return -1;
  }

  if (!geomFileName.size()) {
    std::cerr << "Need geometry." << std::endl;
    return -1;
  }

  if (!mesh.loadobj(geomFileName.c_str())) {
    std::cerr << "Cannot load mesh from:" << geomFileName << std::endl;
    return -2;
  }

  if (bonesFileName.size() && !arm.read(bonesFileName.c_str())) {
    std::cerr << "Error loading bones from: " << bonesFileName << std::endl;
    return -3;
  }

  if (!keys.load(animFileName.c_str())) {
    std::cerr << "Error loading pose keys:" << animFileName << std::endl;
    return -2;
  }

  std::vector<displace_t> disp = dpl::read_dpl(dplFile.c_str(), mesh);
  if (!disp.size()) {
    std::cerr << "Cannot load displacements from:" << dplFile << std::endl;
    return -1;
  }

  nacb::Matrix res = nacb::Matrix::zeros(1, disp.size());
  if (maxNumBasis < 0) {
    maxNumBasis = disp.size();
  }
  maxNumBasis = std::min(maxNumBasis, (int)disp.size());


  // Use the calibration to perform projection and reconstruction.
  if (optind < ac) {
    std::vector<nacb::Matrix> Ps;

    // Load the calibration.
    for (int i=optind; i<ac; i++) 
    {
      nacb::Matrix K(3, 3), E(4, 4), P(3, 4), d(5, 1);
      
      if (!ClbFile::read(av[i], K.data, E.data, P.data, d.data)) {
	std::cerr << "Cannot read calib: " << calibFileName << std::endl;
	return -1;
      }
      P = K * nacb::Matrix::eye(3, 4) * E;
      Ps.push_back(P);

      P.printMatlab("P");
      E.printMatlab("E");
    }


    // Extract vertices that are influenced by enough
    bool extract = false;
    int skip = 4;
    std::vector<int> indices;

    if (extract) {
      std::vector<std::string> muchMotion;
    
      muchMotion.push_back("left_hip");
      muchMotion.push_back("right_hip");
      muchMotion.push_back("left_knee");
      muchMotion.push_back("right_knee");
      muchMotion.push_back("left_toe");
      muchMotion.push_back("right_toe");

      muchMotion.push_back("left_shoulder");
      muchMotion.push_back("right_shoulder");
      muchMotion.push_back("left_elbow");
      muchMotion.push_back("right_elbow");
      muchMotion.push_back("left_wrist");
      muchMotion.push_back("right_wrist");
      muchMotion.push_back("left_hand");
      muchMotion.push_back("right_hand");

      std::set<int> muchMotionIndices;

      for (int i=0; i<muchMotion.size(); i++) {
	muchMotionIndices.insert(arm[muchMotion[i]]->id);
      }
      
      std::vector<int> indices = getAttachedIndices(mesh, muchMotionIndices);
      printf("Attached: %d/%d\n", indices.size(), mesh.vert.size());
    }
    else if (skip > 1) {
      printf("Skipping: %d\n", skip);
      indices = every(mesh, skip);
    }

    if (indices.size()) {
      printf("Extracting mesh part.\n");
      mesh.vert = extractPart(mesh.vert, indices);
      mesh.restVert = extractPart(mesh.restVert, indices);
      disp = extractPart(disp, indices);
    }

    // Zero the displacements
    //disp = std::vector<dpl::displace_t>(disp.size(), dpl::displace_t(mesh.vert.size()));

    // Replace the input with cosine input
    //nacb::Matrix trueB = getCatmullRomBasisMatrix(disp.size(), maxNumBasis + 10);
    //disp = approximateDisplacements(disp, trueB);

    // Get the image projections.
    std::cout << "Getting image projection." << std::endl;
    std::vector<std::vector<std::vector<nacb::Vec2f> > > projections = 
      getImageProjections(Ps, arm, keys, mesh, disp, noiseAmount);

    std::cout << "Got image projection." << std::endl;

    for (int numBasis = maxNumBasis; numBasis >= minNumBasis; numBasis--) {
      nacb::Matrix B(1, 1);

      // Get the appropriate basis.
      if (basisType == "cosine")
	B = getCosineBasisMatrix(disp.size(), numBasis); 
      else if (basisType == "bezier")
	B = getBezierBasisMatrix(disp.size(), numBasis);
      else if (basisType == "catmull")
	B = getCatmullRomBasisMatrix(disp.size(), numBasis);

      srand48(1);
      std::vector<dpl::displace_t> recon = reconstructDisplacements(Ps, arm, keys, mesh, B, projections, plast);
      
      double error = dplDistance(disp, recon);
      std::cout << numBasis << ":" << error << std::endl;
    }
    return 0;
  }

  
  for (int numBasis = maxNumBasis; numBasis >= minNumBasis; numBasis--) {
    nacb::Matrix B(1, 1);

    // Get the appropriate basis.
    if (basisType == "cosine")
      B = getCosineBasisMatrix(disp.size(), numBasis); 
    else if (basisType == "bezier")
      B = getBezierBasisMatrix(disp.size(), numBasis);
    else if (basisType == "catmull")
      B = getCatmullRomBasisMatrix(disp.size(), numBasis);

    PoseKeys currentAnim = keys;
    std::vector<dpl::displace_t> currentDisp = disp;

    // The updating at the bottom tries to fit the animation to any residuals.
    // This does reduce the residual, but only by about a half (across the board).
    for (int k=0; k<numIterations; k++) {
      std::vector<dpl::displace_t> recon = approximateDisplacements(currentDisp, B);
      double error = dplDistance(currentDisp, recon);
      std::cout << numBasis << ":" << error << std::endl;
      res[numBasis - 1] = error;

      dpl::write_dpl("/tmp/recon.dpl", mesh, recon);
      
      if (k + 1 < numIterations) {
	// Get the deformed mesh.
	currentAnim = updateAnimationAndDisplacement(mesh, arm, currentAnim, currentDisp, recon);
	currentAnim.save("/tmp/recon.anim");
      }
    }
  }
  res.printMatlab("res");

  return 0;
}
