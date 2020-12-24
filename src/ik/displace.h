/**
   Common routines for displacements (transform, inverse transform),
   reading, writing, and smoothing.

   These could be merged with some of the routines in gealign subdirectory.
*/
#ifndef IK_DISPLACE_H
#define IK_DISPLACE_H

#include "filters.h"
#include "armature.h"
#include "bone.h"
#include "mesh.h"
#include <vector>
#include <nmath/vec3.h>
#include <nmath/mat4x4.h>

namespace dpl{
  using namespace nacb;
  using namespace std;

  typedef std::vector<Vec3f> displace_t;
  typedef pair<std::vector<Vec3f>, vector<double> > disp_weights_t;

  /**
     write_dpl:  write out a binary form of displacements.
     
     \param fname the output filename
     \param Mesh the mesh for which the displacements belong.
     \param displacements the displacement vector.
  */
  void write_dpl(const char * fname,
		 Mesh & mesh,
		 std::vector<displace_t> & displacements);

  /**
     read_dpl:  write out a binary form of displacements.
     
     \param fname the input filename
     \param Mesh the mesh for which the displacements belong.
     
     \return a vector of displacements.
     \sa write_dpl
  */
  std::vector<displace_t> read_dpl(const char * fname, const Mesh & mesh);

  /**
     extractFromPosed: return the displacements of the mesh, assuming it
     is posed with the given armature.
     
     \return a displacements vector.
   */
  displace_t extractFromPosed(Armature & arm, Mesh & mesh);

  Vec3f transform(std::vector<Mat4x4> & tforms,
		  Mesh::weight_vector_t & weights, 
		  const Vec3f & pt, double w = 0.0);

  Vec3f itransform(std::vector<Mat4x4> & tforms,
		   Mesh::weight_vector_t & weights, 
		   const Vec3f & pt, double w = 0.0);

  void smooth(vector<displace_t> & displacements, int its = 8);

  vector<displace_t> smoothDisplacementsWorld(Armature & arm,
					      PoseKeys & pose,
					      Mesh & mesh,
					      const vector<displace_t> & displacements, 
					      double sigma = 0.1,
					      int nits=2);

};

#endif
