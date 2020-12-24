#ifndef DISPLACE_UV_H
#define DISPLACE_UV_H

#include "ik/mesh.h"
#include <boost/shared_ptr.hpp>
#include <nimage/image.h>
#include <nmath/vec3.h>
#include <nmath/matrix.h>

class BaseGeometry {
 public:
  virtual nacb::Vec3f displace(double x, double y, double z, nacb::Vec3f * normalOut = 0) = 0;
  virtual nacb::Vec3f displaceAndOffset(double x, double y, double z, 
					double ox, double oy, double oz,
					nacb::Vec3f *normalOut = 0,
					nacb::Vec3f *tangentSpaceOut = 0) = 0;
  virtual nacb::Imagef getMask() = 0;
  virtual int getWidth() const = 0;
  virtual int getHeight() const = 0;

  virtual void baseMeshChanged() { }
  
  typedef boost::shared_ptr<BaseGeometry> ptr;
};


class ImageBaseGeometry : public BaseGeometry {
 public:
  nacb::Imagef base;
  nacb::Imagef normals;
  nacb::Imagef mask;

  ImageBaseGeometry(const Mesh& /* mesh, unused */,
		    const nacb::Imagef & _base,
		    const nacb::Imagef & _normals,
		    const nacb::Imagef & _mask,
		    const nacb::Imagef & /* baryImage, unused */,
		    const nacb::Image32& /* _timage, unused */ ) : 
    base(_base), normals(_normals), mask(_mask) { }

  nacb::Vec3f displace(double x, double y, double z, nacb::Vec3f * normalOut = 0){
    nacb::Vec3f co, no;
    base.bilinear(x, y, co.data, 3);
    normals.bilinear(x, y, no.data, 3);

    if (normalOut)
      *normalOut = no;

    return co + no*z;
  }

  /** \brief This is the flowed offset function.  This allows the displacements to
             be represented in a space other than the world coordinates.
   */
  nacb::Vec3f displaceAndOffset(double x, double y, double z, 
				double ox, double oy, double oz,
				nacb::Vec3f *normalOut = 0,
				nacb::Vec3f *tangentSpace = 0) {
    if (tangentSpace) {
      tangentSpace[0] = nacb::Vec3f(1, 0, 0);
      tangentSpace[1] = nacb::Vec3f(0, 1, 0);
      tangentSpace[2] = nacb::Vec3f(0, 0, 1);
    }
    return displace(x, y, z, normalOut) + nacb::Vec3f(ox, oy, oz);
  }

  virtual int getWidth() const {
    return base.w;
  }

  virtual int getHeight() const {
    return base.h;
  }

  nacb::Imagef getMask() {
    return mask;
  }
};


class ImageBaseSkinnedGeometry : public BaseGeometry {
 public:
  ImageBaseSkinnedGeometry(const Mesh& _mesh, 
			   const nacb::Imagef & _base,
			   const nacb::Imagef & _normals,
			   const nacb::Imagef & _mask,
			   const nacb::Imagef & _baryImage,
			   const nacb::Image32 &  _triImage):
    base(_base), normals(_normals), mesh(_mesh), baryImage(_baryImage), triImage(_triImage), mask(_mask) {
    updateTangentSpace();
  }

  nacb::Vec3f displace(double x, double y, double z, nacb::Vec3f * normalOut = 0){
    nacb::Vec3f co(0, 0, 0), no(0, 0, 0);
    nacb::Vec3f bary;

    int tri = triImage((int)round(x), (int)round(y));
    if (tri >= (int)mesh.tris.size())  {
      base.bilinear(x, y, co.data);
      normals.bilinear(x, y, no.data);
      if (normalOut)
	*normalOut = no;
      return co + no*z;
    }
    baryImage.bilinear(x, y, bary.data, 3);

    bary *= (1.0/(bary.x + bary.y + bary.z));

    for (int k=0; k<3; k++) {
      co += mesh.vert[mesh.tris[tri].vi[k]] * bary.data[k];
      no += mesh.norm[mesh.tris[tri].ni[k]] * bary.data[k];
    }
    no.normalize();

    if (normalOut)
      *normalOut = no;

    return co + no*z;
  }

  virtual void baseMeshChanged() { 
    updateTangentSpace();
  }

  void updateTangentSpace() {
    Tu.clear();
    Tv.clear();
    
    Tu = std::vector<nacb::Vec3f>(mesh.vert.size(), nacb::Vec3f(0, 0, 0));
    Tv = std::vector<nacb::Vec3f>(mesh.vert.size(), nacb::Vec3f(0, 0, 0));

    // Compute tangent frame for each triangle.
    for (int ti=0; ti<(int)mesh.tris.size(); ti++) {
      nacb::Vec3f e1 = mesh.vert[mesh.tris[ti].vi[1]] - mesh.vert[mesh.tris[ti].vi[0]];
      nacb::Vec3f e2 = mesh.vert[mesh.tris[ti].vi[2]] - mesh.vert[mesh.tris[ti].vi[0]];
      
      nacb::Vec3f u1 = mesh.tvert[mesh.tris[ti].tci[1]] - mesh.tvert[mesh.tris[ti].tci[0]];
      nacb::Vec3f u2 = mesh.tvert[mesh.tris[ti].tci[2]] - mesh.tvert[mesh.tris[ti].tci[0]];
      
      // The tangent space equations.
      // u1.x * t1.x + u1.y * t1.y = e1;
      // u2.x * t1.x + u2.y * t1.y = e2;

      // Solve for tu, tv
      double det = u1.x*u2.y - u2.x*u1.y;
      nacb::Vec3f tu = u2.y*e1 - u1.y*e2;
      nacb::Vec3f tv = u1.x*e2 - u2.x*e1;

      // Accumulate for each shared vertex.
      for (int k=0; k<3; k++) {
	Tu[mesh.tris[ti].vi[k]] += (tu*(1.0/det));
	Tv[mesh.tris[ti].vi[k]] += (tv*(1.0/det));
      }
    }    

    // Make orthogonal.
    for (int i=0; i<(int)mesh.vert.size(); i++) {
      Tu[i].normalize();
      Tv[i].normalize();

      Tv[i] = mesh.norm[i].cross(Tu[i]);
      Tv[i].normalize();

      Tu[i] = Tv[i].cross(mesh.norm[i]);
      Tu[i].normalize();
    }
  }

  nacb::Vec3f displaceAndOffset(double x, double y, double z, 
				double ox, double oy, double oz,
				nacb::Vec3f *normalOut = 0,
				nacb::Vec3f *tangentSpace = 0) {
    nacb::Vec3f N;
    nacb::Vec3f point = displace(x, y, z, &N);
    nacb::Vec3f U1(1, 0, 0), U2(0, 1, 0);
    if (normalOut)
      *normalOut = N;

    int tri = triImage((int)round(x), (int)round(y));
    nacb::Vec3f bary;

    // Only update TSB if we are inside a triangle.
    if (tri < (int)mesh.tris.size()) {
      baryImage.bilinear(x, y, bary.data, 3);
      bary *= (1.0/(bary.x + bary.y + bary.z));

      U1 = nacb::Vec3f(0, 0, 0);
      U2 = nacb::Vec3f(0, 0, 0);

      for (int k=0; k<3; k++) {
	U1 += (Tu[mesh.tris[tri].vi[k]] * bary.data[k]);
	U2 += (Tv[mesh.tris[tri].vi[k]] * bary.data[k]);
      }
      U1.normalize();
      U2.normalize();
    }

    if (tangentSpace) {
      tangentSpace[0] = U1;
      tangentSpace[1] = U2;
      tangentSpace[2] = N;
    }
    return point + U1*ox + U2*oy + N*oz;
  }  

  virtual int getWidth() const {
    return mask.w;
  }

  virtual int getHeight() const {
    return mask.h;
  }

  nacb::Imagef getMask() {
    return mask;
  }

  nacb::Imagef base;
  nacb::Imagef normals;

  const Mesh& mesh;
  nacb::Imagef baryImage;
  nacb::Image32 triImage;
  nacb::Imagef mask;

  std::vector<nacb::Vec3f> Tu, Tv;
};




#include "ik/mesh.h"
#include <nmath/vec2.h>
#include <nmath/vec3.h>
#include <vector>
#include "DisparitySequence.h"

/**
   This is a warp from uv-coordinates to a particular
   image coordinate frame.  Needs access to the projection
   matrix and also the 3D coordinate mapping from UV to 
   world coordinates.
*/
class DisplaceUV {
 public:
  nacb::Matrix P;
  BaseGeometry::ptr baseGeometry;
  
  DisplaceUV(const BaseGeometry::ptr & _baseGeometry, 
	     const nacb::Matrix & K, const nacb::Matrix & E){
    baseGeometry = _baseGeometry;
    P = K * nacb::Matrix::eye(3, 4) * E;
  }

  inline nacb::Vec3f project(const nacb::Matrix & P, nacb::Vec3f & co, double w = 1.0) const {
    return nacb::Vec3f(P(0, 0)*co.x + P(0, 1)*co.y + P(0, 2)*co.z + P(0, 3)*w,
		       P(1, 0)*co.x + P(1, 1)*co.y + P(1, 2)*co.z + P(1, 3)*w,
		       P(2, 0)*co.x + P(2, 1)*co.y + P(2, 2)*co.z + P(2, 3)*w);
  }

  nacb::Vec2f operator()(float x, float y, float z) const {
    nacb::Vec3f co = baseGeometry->displace(x, y, z);
    nacb::Vec3f proj = project(P, co);
    return nacb::Vec2f(proj.x/proj.z, proj.y/proj.z);
  }

  nacb::Vec2f depthDeriv(float x, float y, float z) const {
    nacb::Vec3f normal(0, 0, 0);
    nacb::Vec3f co = baseGeometry->displace(x, y, z, &normal);
    nacb::Vec3f proj = project(P, co);
    nacb::Vec3f normalProj = project(P, normal, 0.0);

    return nacb::Vec2f(normalProj.x/proj.z - proj.x*normalProj.z/(proj.z*proj.z),
		       normalProj.y/proj.z - proj.y*normalProj.z/(proj.z*proj.z));
  }
};



/**
   This is a warp from uv-coordinates to a particular
   image coordinate frame.  Needs access to the projection
   matrix and also the 3D coordinate mapping from UV to 
   world coordinates.
*/
class DisplaceOffsetUV {
 public:
  nacb::Matrix P;
  BaseGeometry::ptr baseGeometry;
  
  DisplaceOffsetUV(const BaseGeometry::ptr & _baseGeometry, 
		   const nacb::Matrix & K, const nacb::Matrix & E){
    baseGeometry = _baseGeometry;
    P = K * nacb::Matrix::eye(3, 4) * E;
  }

  inline nacb::Vec3f project(const nacb::Matrix & P, const nacb::Vec3f & co, double w = 1.0) const {
    return nacb::Vec3f(P(0, 0)*co.x + P(0, 1)*co.y + P(0, 2)*co.z + P(0, 3)*w,
		       P(1, 0)*co.x + P(1, 1)*co.y + P(1, 2)*co.z + P(1, 3)*w,
		       P(2, 0)*co.x + P(2, 1)*co.y + P(2, 2)*co.z + P(2, 3)*w);
  }

  nacb::Vec2d operator()(float x, float y, float z, const nacb::Vec3d & offs) const {
    nacb::Vec3f co = baseGeometry->displaceAndOffset(x, y, z, offs.x, offs.y, offs.z);
    nacb::Vec3f proj = project(P, co);
    return nacb::Vec2d(proj.x/proj.z, proj.y/proj.z);
  }

  nacb::Vec2d derivatives(float x, float y, float z, const nacb::Vec3d & offs,
			 nacb::Vec2d * oderivs = 0) const {
    nacb::Vec3f tangentSpace[3];
    nacb::Vec3f normal(0, 0, 0);
    nacb::Vec3f co = baseGeometry->displaceAndOffset(x, y, z, offs.x, offs.y, offs.z, &normal, tangentSpace);
    nacb::Vec3f proj = project(P, co);
    nacb::Vec3f normalProj = project(P, normal, 0.0);

    if (oderivs) {
      for (int k=0; k<3; k++) {
	nacb::Vec3f offsProj  = project(P, tangentSpace[k], 0.0);
	oderivs[k] = nacb::Vec2d(offsProj.x/proj.z - proj.x*offsProj.z/(proj.z*proj.z),
				 offsProj.y/proj.z - proj.y*offsProj.z/(proj.z*proj.z));
      }
    }

    return nacb::Vec2d(normalProj.x/proj.z - proj.x*normalProj.z/(proj.z*proj.z),
		       normalProj.y/proj.z - proj.y*normalProj.z/(proj.z*proj.z));
  }
};



// Some utilities for working with displaced meshes.

class DisplacedMesh : public Mesh {
public:
  std::vector<nacb::Vec2<int> > vertexSource;

  DisplacedMesh() {

  }

  void loadSourcesFromTextureCoords(int w, int h) {
    vertexSource = std::vector<nacb::Vec2<int> >(vert.size());
    
    for (int i=0; i<(int)tris.size(); i++) {
      for (int k=0; k<3; k++) {
	vertexSource[tris[i].vi[k]] = nacb::Vec2<int>((int)(tvert[tris[i].tci[k]].x * (w - 1)),
						      (int)(tvert[tris[i].tci[k]].y * (h - 1)));
      }
    }
  }

  void displace(const BaseGeometry::ptr & baseGeometry,
		nacb::Imagef & depth) {
    for (int i=0; i<(int)vert.size(); i++) {
      const nacb::Vec2<int> & co2 = vertexSource[i];
      vert[i] = baseGeometry->displace(co2.x, co2.y, depth(co2.x, co2.y));
    }
    initNormals();
  }

  void displaceAndOffset(const BaseGeometry::ptr & baseGeometry,
			 const nacb::Imagef & depth) {
    for (int i=0; i<(int)vert.size(); i++) {
      const nacb::Vec2<int> & co2 = vertexSource[i];
      vert[i] = baseGeometry->displaceAndOffset(co2.x, co2.y, depth(co2.x, co2.y),
						depth(co2.x, co2.y, 1), 
						depth(co2.x, co2.y, 2), 
						depth(co2.x, co2.y, 3));

      if (!((isnormal(vert[i].x) || vert[i].x == 0) && 
	    (isnormal(vert[i].y) || vert[i].y == 0) && 
	    (isnormal(vert[i].z) || vert[i].z == 0))) {
	vert[i] = restVert[i];
	printf("not normal: %f %f %f\n", vert[i].x, vert[i].y, vert[i].z);
      }
    }
    initNormals();
  }
};


// Assumes a GL context exists.
template <class BaseGeometryType>
BaseGeometry::ptr createBaseGeometryFromMesh(Mesh & mesh, int w, int h, 
					     DisplacedMesh & submesh);


nacb::Imagef getWeightImage(const DisparitySequence & seq,
			    BaseGeometry::ptr & baseGeometry, 
			    DisplacedMesh * submesh = 0,
			    bool occlusion = false);

#endif // DISPLACE_UV_H
