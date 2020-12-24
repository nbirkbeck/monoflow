#include "GroundTruthLib.h"
#include <nmath/matrix.h>

#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"
#include "autorecon/stereo/simplesequence.h"

using namespace nacb;


std::vector<Triangle> getTriangles(const Mesh & mesh){
  std::vector<Triangle> tris(mesh.tris.size());

  for(int ti=0; ti<(int)tris.size(); ti++){
    tris[ti] = Triangle(mesh.vert[mesh.tris[ti].vi[0]],
			mesh.vert[mesh.tris[ti].vi[1]],
			mesh.vert[mesh.tris[ti].vi[2]]);
  }

  return tris;
}


double getClosestTriangle(const std::vector<Triangle> & triangles,
			  const Vec3f & cc,
			  const Vec3f & dir,
			  int & which,
			  Vec3f & weights){
  double depth = -1;

  which = -1;

  for(int i=0; i<(int)triangles.size(); i++){
    double dray = 0;
    double a, b, c;

    if(triangles[i].intersect(cc, dir, &dray, &a, &b, &c)){
      if((depth < 0) || (dray < depth)){
	depth = dray;
	weights = Vec3f(a, b, c);
	which = i;
      }
    }
  }
  
  return depth;
}


Vec3f getBaryVertex(const Triangle & triangle, const Vec3f & w){
  return triangle.v1*w.x  + triangle.v2*w.y + triangle.v3*w.z;
}


void averageOutside(const std::vector<Triangle> & triangles,
		    const nacb::Image<uint32_t> & itri,
		    nacb::Imagef & depth,
		    const Vec3f & average,
		    int maxits) {
  // Compute average in border (and propagate it)
  nacb::Image<uint32_t> itriUse = itri.copy();

  for (int its=0; its<maxits; its++) {
    nacb::Image<uint32_t> itriNew = itriUse;

    nacb::Imagef::iterator last = depth.end();
    
    for(nacb::Imagef::iterator it = depth.begin(); it != last; ++it){
      if(itriUse(it.x, it.y) == triangles.size()) {
	int neigh[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	Vec3f lavg(0, 0, 0);
	int navg = 0;

	for (int k=0; k<4; k++){
	  int x2 = it.x + neigh[k][0];
	  int y2 = it.y + neigh[k][1];
	  
	  if (x2>=0 && y2>=0 && x2<depth.w && y2<depth.h &&
	      itriUse(x2, y2) < triangles.size()) {
	    for(int k=0; k<depth.nchannels; k++)
	      lavg.data[k] += depth(x2, y2, k);
	    navg++;
	  }
	}

	if (navg > 0) {
	  for(int k=0; k<depth.nchannels; k++)
	    depth(it.x, it.y, k) = lavg.data[k]/navg;

	  itriNew(it.x, it.y) = 0;
	}
      }
    }
    
    itriUse = itriNew.copy();
  }
  
  nacb::Imagef::iterator last = depth.end();

  for(nacb::Imagef::iterator it = depth.begin(); it != last; ++it){
    if(itriUse(it.x, it.y) == triangles.size()) {
      for (int k=0; k<depth.nchannels; k++)
	depth(it.x, it.y, k) = average.data[k];
    }
  }
}

Imagef computeImageFlow(const nacb::Image8 & image,
			const std::vector<Triangle> & tri1,
			const std::vector<Triangle> & tri2,
			const nacb::Image<uint32_t> & itri,
			const nacb::Imagef & weights){
  nacb::Imagef flow(image.w, image.h, 3);

  flow = 0;
  for(int y=0; y<image.h; y++){
    for(int x=0; x<image.w; x++){
      if(itri(x, y) < tri1.size()){
	Vec3f bary = Vec3f(weights(x, y, 0), weights(x, y, 1), weights(x, y, 2));
	Vec3f p2 = getBaryVertex(tri2[itri(x, y)], bary);
	Vec3f p1 = getBaryVertex(tri1[itri(x, y)], bary);
	
	flow(x, y, 0) = p2.x - p1.x;
	flow(x, y, 1) = p2.y - p1.y;
	flow(x, y, 2) = p2.z - p1.z;
      }
    }
  }

  averageOutside(tri1, itri, flow, Vec3f(0, 0, 0), 10);

  return flow;
}


Imagef computeImageDisparity(nacb::Image8 & image, Matrix & P, 
			     const std::vector<Triangle> & triangles,
			     nacb::Image<uint32_t> & itri,
			     nacb::Imagef & weights){
  nacb::Imagef depth(image.w, image.h, 1);
  
  Matrix KRinv = P.submatrix(0, 0, 3, 3).inverse();
  Matrix Kt = P.submatrix(0, 3, 3, 1);

  Vec3f cc = backProject(KRinv, Kt, 0.f, 0.f, 0.f);
  nacb::Image8::iterator last = image.end();

  itri = nacb::Image<uint32_t>(image.w, image.h, 1);
  weights = nacb::Imagef(image.w, image.h, 3);

  double avg = 0;
  int navg = 0;

  for(nacb::Image8::iterator it = image.begin(); it != last; ++it){
    Vec3f ray = backProject(KRinv, Kt, (float)it.x, (float)it.y, 1.0f) - cc;
    double zscale = ray.normalize();
    bool good = false;

    if((*it)[3] >= 120){
      int tri;
      Vec3f wts;
      double dray = getClosestTriangle(triangles, cc, ray, tri, wts)/zscale;

      if(dray > 0){
	depth(it.x, it.y) = 1.0/dray;

	//float sum = wts.x + wts.y + wts.z;
	//Vec3f v = getBaryVertex(triangles[tri], wts) - cc;
	//v.normalize();
	//printf("%f\n", v.dot(ray)); 

	weights(it.x, it.y, 0) = wts.x;
	weights(it.x, it.y, 1) = wts.y;
	weights(it.x, it.y, 2) = wts.z;

	itri(it.x, it.y) = tri;

	good = true;

	avg += depth(it.x, it.y);
	navg ++;
      }
    }
    if(!good){
      depth(it.x, it.y) = 0.0f;
      itri(it.x, it.y) = triangles.size();
    }
  }

  // Set background to average...this was a bad idea.  Need to set the depth to the
  avg /= navg;

  averageOutside(triangles, itri, depth, Vec3f(avg, avg, avg), 10);
  
  return depth;
}

