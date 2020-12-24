/** Common routines for dealing with problems involving a displaced offset proxy
    that uses a basis.

    Neil Birkbeck,
*/

#ifndef DISPLACE_OFFSET_BASIS
#define DISPLACE_OFFSET_BASIS 1

#include "DisplaceUV.h"
#include "BaseMeshAnimation.h"
#include "FlowBasis3D.h"
#include <nimage/image.h>
#include <vector>

namespace DisplaceOffsetBasis {
  
  inline nacb::Imagef warpBack(const std::vector<nacb::Imagef> & images, 
			       const FlowBasis3D::ptr & basis,
			       BaseMeshAnimation& animation,
			       DisplaceOffsetUV & warp,
			       const nacb::Imagef & D,
			       int t = 0){
    nacb::Imagef warped(D.w, D.h, images[t].nchannels);
    warped = 0;

    animation.setTime(t);

    for (int y=0; y<D.h; y++) {
      for (int x=0; x<D.w; x++) {
	double displace;
	Vec3f offs = (*basis)(displace, x, y, D, t);
	Vec2f co = warp(x, y, displace, offs);
      
	for (int k=0; k<images[t].nchannels; k++) {
	  warped(x, y, k) = images[t].bilinear(co.x, co.y, k);
	}
      }
    }
    return warped;
  }


  inline std::vector<std::pair<int, int> > buildStereoPairs(int n) {
    std::vector<std::pair<int, int> > pairs;
    
    for (int i=0; i<n; i++) {
      for (int j=i+1; j<n; j++) {
	pairs.push_back(make_pair(i, j));
      }
    }
    return pairs;
  }


  inline void initImagesWarpsAndWeights(std::vector<DisparitySequence>& seqs,
				 FlowBasis3D::ptr& basis,
				 BaseMeshAnimation& meshAnimation,
				 DisplacedMesh& dispMesh,
				 BaseGeometry::ptr& baseGeometrySubdiv,
				 const std::vector<int>& times,
				 nacb::Imagef& depthOffset,
				 std::vector<std::vector<nacb::Imagef> >&  images,
				 std::vector<DisplaceOffsetUV>&     warps,
				 std::vector<std::vector<nacb::Imagef> >&   weights, 
				 bool saveWeights = false) 
  {
    // Initialize the warps and the image list.
    images.clear();
    warps.clear();
    weights.clear();
    
    for (int i=0; i < (int)seqs.size(); i++) {
      int w = seqs[i].image.w, h = seqs[i].image.h;
      vector<Imagef> imagesi;
      vector<Imagef> weightsi;
      
      for (int ti = 0; ti < (int)times.size(); ti++) {
	int t1 = times[ti];
	
	meshAnimation.setTime(t1);
	dispMesh.displaceAndOffset(baseGeometrySubdiv, basis->getDisplacementMap(depthOffset, t1));
	
	seqs[i].load(t1);
	seqs[i].image = seqs[i].image.resize(w, h);
	imagesi.push_back(seqs[i].image);
	
	weightsi.push_back(getWeightImage(seqs[i], baseGeometrySubdiv, &dispMesh, true));
	if (saveWeights)
	  weightsi.back().write((boost::format("/tmp/weights-%04d-%02d.png") % i % t1).str().c_str());
      }
      
      warps.push_back(DisplaceOffsetUV(baseGeometrySubdiv, seqs[i].A, seqs[i].E));      
      images.push_back(imagesi);
      weights.push_back(weightsi);
    }
  }

  //!\brief Normalize the weight images so that all pixels have roughly the same strength (help 
  //        equalize the regularization over the surface)
  // Didn't seem to help all that much.
  inline void normalizeWeights(std::vector<nacb::Imagef>&   weights) {
    nacb::Imagef sum(weights[0].w, weights[0].h, 1);
    // nacb::Imagef maximum(weights[0].w, weights[0].h, 1);
    sum = 0.1; // Use an epsilon to avoid divide by zero.
    // maximum = 0.1;

    for (int j=0; j < (int)weights.size(); j++) {
      sum += weights[j];

      /*
      for (int y=0; y<sum.h; y++)
	for (int x=0; x<sum.w; x++) {
	  maximum(x, y) = max(maximum(x, y), weights[j](x, y));
	}
      */
    }

    for (int j=0; j < (int)weights.size(); j++) {
      for (int y=0; y < sum.h; y++)
	for (int x=0; x < sum.w; x++) {
	  weights[j](x, y) /= sum(x, y);
	  weights[j](x, y) *= weights.size();
	}
    }
  }
};

#endif // DISPLACE_OFFSET_BASIS
