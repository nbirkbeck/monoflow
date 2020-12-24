/**
   Some helper functions to generate ground truth from a rendered sequence.

   Neil Birkbeck, 2009-2011
 */
#ifndef GROUND_TRUTH_LIB_H
#define GROUND_TRUTH_LIB_H

#include "ik/triangle.h"
#include "ik/mesh.h"
#include <vector>
#include <nmath/vec3.h>
#include <nmath/matrix.h>
#include <nimage/image.h>

std::vector<Triangle> getTriangles(const Mesh & mesh);
double getClosestTriangle(const std::vector<Triangle> & triangles,
			  const nacb::Vec3f & cc,
			  const nacb::Vec3f & dir,
			  int & which,
			  nacb::Vec3f & weights);

nacb::Vec3f getBaryVertex(const Triangle & triangle, const nacb::Vec3f & w);

void averageOutside(const std::vector<Triangle> & triangles,
		    const nacb::Image<uint32_t> & itri,
		    nacb::Imagef & depth,
		    const nacb::Vec3f & average,
		    int maxits);

nacb::Imagef computeImageFlow(const nacb::Image8 & image,
			      const std::vector<Triangle> & tri1,
			      const std::vector<Triangle> & tri2,
			      const nacb::Image<uint32_t> & itri,
			      const nacb::Imagef & weights);

nacb::Imagef computeImageDisparity(nacb::Image8 & image, nacb::Matrix & P, 
				   const std::vector<Triangle> & triangles,
				   nacb::Image<uint32_t> & itri,
				   nacb::Imagef & weights);


#endif
