/**
   Use the long-range tracks (in 2D) to fit the armature motion.

   The long range tracks are just a set of correspondences over
   a few frames.
   
   Assuming the mesh in the first frame is accurate, we want
   the motion over the sequence (from one frame to the next)
   to obey the point constraints.

   E.g., for time 1 to time 2:

   \sum_k |(I(P_i(p(\theta_2))) - I(P_i(p(\theta_1)))) - [u,v]|^2

   Where p is a bary centric point on the mesh
   E.g., just a non-linear equation system in the joint angles.

   Solve for angles (given the theta_1)

   Warning: the analytical derivative is wrong.

   Jan 2011
*/

#include "DisparitySequence.h"
#include "GroundTruthLib.h"
#include "autorecon/recon_globals.h"
#include "autorecon/apps/manual/observedpoint.h"
#include <nmisc/commandline.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include "ik/mesh.h"
#include "ik/armature.h"
#include "ik/triangle.h"
#include <set>

/** \!brief  A bary centric point on a triangle.
*/
struct BaryPoint {

  nacb::Vec3f getValue(const Mesh& mesh) const {
    nacb::Vec3f co(0, 0, 0);
    for (int k=0; k<3; k++) {
      co += mesh.vert[mesh.tris[triangle].vi[k]] * weight.data[k];
    }
    return co;
  }

  nacb::Vec3f weight; //!< The bary weight.
  int       triangle; //!< The triangle.
};



/** \!brief A correspondence between a mesh point and an image offset.
*/
struct Correspondence {
  Correspondence(const BaryPoint& p = BaryPoint(),
		 const nacb::Vec2f& co = nacb::Vec2f(0, 0),
		 int i = 0): 
    point(p), coordinate(co), image(i) { }

  BaryPoint         point; //!< The point
  nacb::Vec2f  coordinate; //!< Image coordinate.
  int image; //!< image index.
};



/**
   The context of the problem.   Holds the mesh and the cameras.
*/
class TemporalTrackProblem {
public:
  TemporalTrackProblem() { }

  // Evaluate objective function.
  double objective(int m, int n, const double* x, double* fvec) {
    assert(m == correspondences.size()*2);

    arm.setParams(x);
    arm.animate(mesh);

    double residual = 0.0;

    // Each correspondence has two equations.
    for (int i=0; i<correspondences.size(); i++) {
      Correspondence& co = correspondences[i];
      nacb::Vec2f p = project(Ps[co.image], co.point.getValue(mesh));

      fvec[2*i    ] = p.x - co.coordinate.x;
      fvec[2*i + 1] = p.y - co.coordinate.y;

      residual += fvec[2*i]*fvec[2*i];
      residual += fvec[2*i + 1]*fvec[2*i + 1];
    }
    return sqrt(residual / correspondences.size());
  }

  // Compute the jacobian (Warning: this derivative is wrong); use the numerical derivative.
  void jacobian(int m, int n, const double* x, const double* fvec, double * fjac, int ldfjac) {
    assert(m == correspondences.size()*2);

    arm.setParams(fvec);
    arm.animate(mesh);

    double * J = arm.meshJacobian(mesh, fvec, n, false, &jcache);

    memset(fjac, 0, sizeof(double)*ldfjac*n);

    // Each correspondence has two equations.
    for (int i=0; i<correspondences.size(); i++) {
      Correspondence& co = correspondences[i];

      float z;
      nacb::Vec3f vert = co.point.getValue(mesh);
      nacb::Vec2f p = project(Ps[co.image], vert, &z);

      nacb::Vec3f p3(p.x*z, p.y*z, z);

      // Each correspondence is a function of three bary points.
      for (int k=0; k<3; k++) {
	double weight = co.point.weight.data[k];
	int point = mesh.tris[co.point.triangle].vi[k];

	// Each of the bary points is possibly a function of all joint angles.
	for (int j=0; j<n; j++) {
	  // Grab the derivative of the point w.r.t to the angle.
	  double * dPdj = J + j*3*mesh.vert.size() + 3*point;

	  // This was wrong before because it should be multiplying the dPI/dx by P before
	  nacb::Vec3d P_dPdj(0, 0, 0);
	  for (int ii=0; ii<3; ii++) {
	    for (int jj=0; jj<3; jj++) {
	      P_dPdj[ii] += Ps[co.image](ii, jj) * dPdj[jj];
	    }
	  }
	  fjac[ldfjac*j + 2*i    ] += weight*(P_dPdj[0]/z - p3.x*P_dPdj[2]/(z*z));
	  fjac[ldfjac*j + 2*i + 1] += weight*(P_dPdj[1]/z - p3.y*P_dPdj[2]/(z*z));
	}
      }
    }

    delete [] J;
  }

  int getNumFunctions() const {
    return correspondences.size()*2;
  }

  std::set<int> getTrianglesWithoutTracks() {
    std::set<int> withTrack;

    for (int i=0; i<correspondences.size(); ++i) {
      withTrack.insert(correspondences[i].point.triangle);
    }

    std::set<int> result;
    for (int i=0; i<mesh.tris.size(); i++) {
      if (!withTrack.count(i))
	result.insert(i);
    }
    return result;
  }

  //!\brief Solve for the joint angles given the current settings.
  nacb::Matrix solve() {
    nacb::Matrix soln = initAngles.copy();
    nacb::Matrix fvec(getNumFunctions(), 1);
    
    double before = objective(fvec.m, soln.n, soln.data, fvec.data);
    //nacb::Matrix::lsqnonlin(s_objective, (void*)this, soln, fvec);
    nacb::Matrix::lsqnonlin(s_objective2, (void*)this, soln, fvec);
    double after =  objective(fvec.m, soln.n, soln.data, fvec.data);
    printf("res: %f, %f\n", before, after);
    return soln;
  }
  
  //!\brief  The ugly objective (used in the lsqnonlin api).
  static int s_objective2(void * p, int m, int n, 
			 const double * x, double * fvec,
			 int iflag) {
    ((TemporalTrackProblem*)p)->objective(m, n, x, fvec);
    return 0;
  }


  

  //!\brief  The ugly objective (used in the lsqnonlin api).
  static int s_objective(void * p, int m, int n, 
			 const double * x, double * fvec,
			 double * fjac,
			 int ldfjac,
			 int iflag) {
    if (iflag == 1) {
      ((TemporalTrackProblem*)p)->objective(m, n, x, fvec);
    }
    else if (iflag == 2) {
      ((TemporalTrackProblem*)p)->jacobian(m, n, x, fvec, fjac, ldfjac);
    }
    return iflag;
  }

  //!\brief Setup the internal angles to be the defaults.
  void useDefaultAngles() {
    initAngles = nacb::Matrix(arm.getNParams(), 1);
    arm.getParams(initAngles.data);
  }

  //!\brief Clear the correspondenecs.
  void clearCorrespondences() {
    correspondences.clear();
  }
  
  //!\brief Use the cameras from the given sequences.
  void setCameras(const std::vector<DisparitySequence>& seqs) {
    Ps.clear();

    for (int i=0; i<seqs.size(); i++) {
      Ps.push_back(seqs[i].P);
    }
  }

  void animateMesh() {
    arm.animate(mesh);
  }

public:
  JacobianCache jcache;

  Armature arm; //!< The armature.
  Mesh    mesh; //!< A mesh.

  std::vector<Correspondence> correspondences; //!< Point correspondences.
  std::vector<nacb::Matrix>                Ps; //!< Camera matrices.

  nacb::Matrix initAngles;
};



int main(int ac, char * av[]) {
  nacb::CommandLine cline;
  std::string geomFile;
  std::string bonesFileName;
  std::string pointsFile;

  int numFrames = 10;

  cline.registerOption("geom", "The geometry file", &geomFile, 'g');
  cline.registerOption("bones", "Bones file name", &bonesFileName, 0);
  cline.registerOption("numFrames", "The number of frames to use.", &numFrames);
  cline.registerOption("points", "The file with the points.", &pointsFile);
  cline.parse(ac, av);


  if (optind == ac) {
    std::cout << "Need to provide image sequences." << std::endl;
    return -1;
  }

  std::vector<DisparitySequence> seqs = 
    DisparitySequence::loadSequences(av + optind, ac - optind);

  std::vector<ObservedPoints> points(seqs.size());

  for (int i=0; i<seqs.size(); i++) {
    if (!points[i].load((boost::format(pointsFile) % i).str().c_str())) {
      std::cout << "Cannot read points file: " << pointsFile << " " << i << std::endl;
      return -1;
    }
  }

  TemporalTrackProblem problem;

  if (!problem.mesh.loadobj(geomFile.c_str())) {
    std::cout << "Cannot read the geometry:" << geomFile << std::endl;
    return -1;
  }

  if (!problem.arm.read(bonesFileName.c_str())) {
    std::cout << "Cannot read the bones:" << bonesFileName << std::endl;
    return -1;
  }

  problem.setCameras(seqs);
  problem.arm.useEuclideanRoot();
  problem.useDefaultAngles();


  PoseKeys keys(problem.arm.getBoneNames());
  keys.append(&problem.arm);

  // For each frame, for each image, find correspondences, solve for new angles.
  for (int frame=0; frame < (numFrames-1); frame++) {
    std::vector<Triangle> triangles = getTriangles(problem.mesh);

    // Reset the problem constraints.
    problem.animateMesh();
    problem.clearCorrespondences();

    // For each image sequence.
    for (int si=0; si<seqs.size(); si++) {
      std::cout << "Getting correspondences for:" << si << std::endl;
      std::vector<ObservedPoint> ops = points[si];

      nacb::Vec3f cc = seqs[si].cc;
      int numCandidates = 0;
      int numFound = 0;
      double averageMotion = 0.0;

      for (int pointIndex=0; pointIndex<ops.size(); pointIndex++) {
	// If correspondence exists frame and (frame+1)
	if (ops[pointIndex].observedBy(frame) &&
	    ops[pointIndex].observedBy(frame + 1)) {
	  nacb::Vec2f p1 = ops[pointIndex].getObservation(frame);
	  nacb::Vec2f p2 = ops[pointIndex].getObservation(frame + 1);

	  // Get back-projected ray
	  nacb::Vec3f dir = backProject(seqs[si].P, p1.x, p1.y, 1.f);
	  dir = dir - cc;
	  dir.normalize();

	  
	  // Find bary-centric coordinates on most recent mesh.
	  BaryPoint point;
	  if (getClosestTriangle(triangles, cc, dir, point.triangle, point.weight) > 0) {
	    numFound++;
	    problem.correspondences.push_back(Correspondence(point, p2, si));
	    averageMotion += (p1 - p2).len();
	  }
	  numCandidates++;
	}
      }
      std::cout << "Found:" << numFound << " Total candidates:" << numCandidates << std::endl;
      std::cout << "Average motion:" << averageMotion/numCandidates << std::endl;
    }
    std::cout << "Total correspondences:" << problem.correspondences.size() << std::endl;

    std::set<int> noTracks = problem.getTrianglesWithoutTracks();
    std::set<int>::iterator it(noTracks.begin());
    const int numCons = 100;

    std::cout << "Missing tracks:" << noTracks.size() << " these will be forced to remain stationary." << std::endl;

    for ( ; it != noTracks.end(); ++it) {
      for (int si=0; si<seqs.size(); si++) {
	for (int j=0; j<numCons; j++) {
	  BaryPoint point;
	  point.triangle = *it;
	  point.weight = nacb::Vec3f(0.333f, 0.333f, 0.333f);

	  nacb::Vec3f P = point.getValue(problem.mesh);
	  nacb::Vec2f p = project(seqs[si].P, P);

	  problem.correspondences.push_back(Correspondence(point, p, si));
	}
      }
    }

    // Solve for angles.
    nacb::Matrix angles = problem.solve();

    // Add the latest angles to the set of pose keys.
    problem.arm.setParams(angles.data);
    keys.append(&problem.arm);
    problem.initAngles = angles;
  }

  keys.save("/tmp/long-range.anim");
  return 0;
}
