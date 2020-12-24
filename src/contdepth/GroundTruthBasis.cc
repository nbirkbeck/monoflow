/*
  Compute the ground truth for the basis flow.

  Neil Birkbeck, Jan 2011
*/
#include <algorithm>
#include <vector>

#include "ik/mesh.h"
#include "ik/triangle.h"

#include "DisplaceUV.h"
#include "GroundTruthLib.h"

#include "utigl/glwindow.h"
#include <nmisc/commandline.h>
#include "TimeVaryingDisplacedMesh.h"

using namespace nacb;


// Compute the closest triangle for each point in the base mesh.
Imagef computeMeshDisplacement(const BaseGeometry::ptr& baseGeom,
			       const std::vector<Triangle> & triangles,
			       nacb::Image<uint32_t> & itri,
			       nacb::Imagef & weights) {
  double avg = 0;
  int navg = 0;
  nacb::Imagef depth(baseGeom->getWidth(), baseGeom->getHeight(), 1);
  depth = 0;

  weights = nacb::Imagef(depth.w, depth.h, 3);
  weights = 0;

  itri = nacb::Image<uint32_t>(depth.w, depth.h, 1);
  nacb::Imagef mask = baseGeom->getMask();

  
  for (int y=0; y<depth.h; y++) {
    for (int x=0; x<depth.w; x++) {
      if (mask(x, y) > 0) {
	Vec3f ray;
	Vec3f cc = baseGeom->displace(x, y, 0, &ray);
	bool good = false;
	
	Vec3f wts;
	double dRay = 1e10;
	int tri = -1;
	
	for (int ti=0; ti<(int)triangles.size(); ti++) {
	  double a, b, c;
	  double d;
	  if (triangles[ti].intersect(cc, ray, &d, &a, &b, &c)){
	    if (fabs(d) < fabs(dRay)) {
	      dRay = d;
	      wts.x = a;
	      wts.y = b;
	      wts.z = c;
	      tri = ti;
	    }
	  }
	}

	if (tri >= 0) {
	  depth(x, y) = dRay;
	  
	  //float sum = wts.x + wts.y + wts.z;

	  weights(x, y, 0) = wts.x;
	  weights(x, y, 1) = wts.y;
	  weights(x, y, 2) = wts.z;
	  
	  itri(x, y) = tri;
	  
	  good = true;
	  
	  avg += depth(x, y);
	  navg ++;
	}
      
	if(!good){
	  itri(x, y) = triangles.size();
	}
      }
      else
	itri(x, y) = triangles.size();
    }
  }
  avg /= navg;

  averageOutside(triangles, itri, depth, avg, 15);
  return depth;
}



int main(int ac, char * av[]){
  nacb::CommandLine cline;
  std::string geomFile, animFile, bonesFile;
  std::string truthFile;
  std::string basisType;
  int referenceTimeIndex = 0;
  int index1 = 0;
  int index2 = 4;
  int numBasis = 4;
  int texw = 512, texh = 512;
  int depthOnly = 0;

  cline.registerOption("geom", "Geometry base with %%d in it", &geomFile, 'g');
  cline.registerOption("bones", "The bones file", &bonesFile, 'b');
  cline.registerOption("anim", "Geometry base with %%d in it", &animFile, 'a');
  cline.registerOption("truth", "The ground truth geometry with %%d in it", &truthFile, 't');
  cline.registerOption("basisType", "The basis type.", &basisType, 0);
  cline.registerOption("index1", "The index of the first frame.", &index1, 0);
  cline.registerOption("index2", "The index of the last time frame.", &index2, 0);
  cline.registerOption("numBasis", "The number of basis elements.", &numBasis, 0);
  cline.registerOption("texw", "The texture width.", &texw, 0);
  cline.registerOption("texh", "The texture height.", &texh, 0);
  cline.registerOption("depthOnly", "Only reconstruct depth (not flow).", &depthOnly, 0);
  cline.parse(ac, av);

  // Need a window to create base geometry.
  GLWindow win(1024, 1024);

  TimeVaryingDisplacedMesh timeVaryingMesh;

  if (!timeVaryingMesh.loadMesh(geomFile)) {
    std::cerr << "Cannot load base mesh." << std::endl;
    return -1;
  }

  if (!timeVaryingMesh.loadKinematics(bonesFile, animFile)) {
    std::cout << "Warning: cannot load bones or animation.\n";
  }

  FlowBasis3D::ptr basis(new StandardFlowBasis3D(referenceTimeIndex));
  if (basisType == "cosine") 
    basis = FlowBasis3D::ptr(new CosineFlowBasis3D(referenceTimeIndex, numBasis, (index2 - index1 + 1)));
  else if (basisType == "perframe") {
    basis = FlowBasis3D::ptr(new PerFrameFlowBasis3D(referenceTimeIndex, (index2 - index1 + 1)));
  }
  else if (basisType == "cosinedepth") {
    basis = FlowBasis3D::ptr(new CosineDepthOffset(referenceTimeIndex, numBasis, (index2 - index1 + 1)));
  }
  if (!timeVaryingMesh.setBasis(basis)) {
    std::cerr << "Cannot load basis file.\n";
    return -1;
  }
  // Resize the basis data.
  timeVaryingMesh.resizeBasisData(texw, texh);

  Mesh groundTruthMesh;
  groundTruthMesh.loadobj((boost::format(truthFile) % referenceTimeIndex).str().c_str());

  std::vector<Triangle> triangles = getTriangles(groundTruthMesh);
  nacb::Imagef weights;
  nacb::Image<uint32_t> intersects;
  nacb::Imagef depth = computeMeshDisplacement(timeVaryingMesh.getBaseGeometry(), triangles, intersects, weights);
  
  depth.getNormalizedImage().write("/tmp/_gt.png");


  // Compute the correspondences for each mesh.
  std::vector<nacb::Matrix> correspondences;
  nacb::Image8 mask(intersects.w, intersects.h, 1);
  mask = 0;

  bool hasDisp = basis->hasDisplacementComponent();

  for (int y=0; y<intersects.h; y++) {
    for (int x=0; x<intersects.w; x++) {
      if (intersects(x, y) < triangles.size()) {
	nacb::Matrix M = nacb::Matrix::zeros(4*(index2 - index1 + 1), 1);

	// Set the depth observation (over time)
	if (hasDisp) {
	  for (int t=0; t<M.m; t+=4)
	    M[t] = depth(x, y);
	}
	correspondences.push_back(M);
	mask(x, y) = 255;
      }
    }
  }
  mask.save("/tmp/mask.png");


  for (int t=index1; t<=index2; t++) {
    Mesh displacedTruth;
    if (!displacedTruth.loadobj((boost::format(truthFile) % t).str().c_str())) {
      std::cout << "Cannot load mesh from time: " << t << std::endl;
      return 0;
    }

    // Set the pose of the mesh, so the displacements are relative to animation.
    timeVaryingMesh.setTimePoseOnly(t);
    
    nacb::Image<uint32_t> intersectsThisTime;
    nacb::Imagef depthThisTime;

    // If we are just reconstructing depth maps, we will have to find the displacement map again:
    if (depthOnly) {
      std::vector<Triangle> triangles = getTriangles(displacedTruth);
      nacb::Imagef weightsThisTime;
      depthThisTime = computeMeshDisplacement(timeVaryingMesh.getBaseGeometry(), triangles, intersectsThisTime, weightsThisTime);
    }


    // For each pixel that is on, compute the displacement of the vertex.
    int mi = 0;
    nacb::Matrix M(3, 3);

    for (int y=0; y<intersects.h; y++) {
      for (int x=0; x<intersects.w; x++) {
	nacb::Vec3f vert(0, 0, 0);
	//nacb::Vec3f restVert(0, 0, 0);

	if (intersects(x, y) < triangles.size()) {
	  if (depthOnly) {
	    vert = timeVaryingMesh.getBaseGeometry()->displace(x, y, depthThisTime(x, y));
	  }
	  else {
	    int ti = intersects(x, y);

	    for (int k=0; k<3; k++) {
	      vert += displacedTruth.vert[displacedTruth.tris[ti].vi[k]] * weights(x, y, k);
	      //restVert += groundTruthMesh.vert[groundTruthMesh.tris[ti].vi[k]] * weights(x, y, k);
	    }
	  }

	  // The actual rest vert should come from the baseGeometry
	  nacb::Vec3f normal;
	  nacb::Vec3f tangentSpace[3];

	  // If there is no displacement component, allow the offsets to model it (e.g., use base mesh as rest vert)
	  nacb::Vec3f restVert = timeVaryingMesh.getBaseGeometry()->displaceAndOffset(x, y, hasDisp ? depth(x, y): 0.0,
										      0.f, 0.f, 0.f,
										      &normal,
										      tangentSpace);	  
	  // This is now the offset (in world coordinates).
	  vert -= restVert;
	  
	  nacb::Matrix M(3, 3);
	  for (int j=0; j<3; j++) {
	    M(0, j) = tangentSpace[j].x;
	    M(1, j) = tangentSpace[j].y;
	    M(2, j) = tangentSpace[j].z;
	  }

	  // Put it into tangent space
	  nacb::Matrix b = vert;
	  vert = nacb::Matrix::LlinLeastSq(M, b);

	  //int i = t - index1;
	  correspondences[mi][t*4 + 1] = vert.data[0];
	  correspondences[mi][t*4 + 2] = vert.data[1];
	  correspondences[mi][t*4 + 3] = vert.data[2];
	  mi++;
	}
      }
    }
  }


  // For pixel that has a correspondence, now compute the basis coefficients. 
  int mi = 0;
  std::vector<int> times;
  for (int t=index1; t<=index2; t++) 
    times.push_back(t);

  nacb::Matrix B = basis->getBasisMatrix(times).transpose();
  nacb::Imagef basisRep(depth.w, depth.h, B.n);
  basisRep = 0;

  int numPoints = 0;
  double pointRes = 0;
  double reconRes = 0;

  for (int y=0; y<intersects.h; y++) {
    for (int x=0; x<intersects.w; x++) {
      if (intersects(x, y) < triangles.size()) {
	nacb::Matrix coeff = nacb::Matrix::linLeastSq(B, correspondences[mi]);
	nacb::Matrix res = B*coeff - correspondences[mi];

	for (int k=0; k<basisRep.nchannels; k++) {
	  basisRep(x, y, k) = coeff[k];
	}
	pointRes += sqrt(res.dot(res) / times.size());
	numPoints++;

	// This other residual provides a way to test that the basismatrix is the same
	// as what is given by operator().
	double recon = 0.0;
	for (int t=index1; t<=index2; t++) {
	  double displace;
	  Vec3f offs = (*basis)(displace, x, y, basisRep, t);
	  int ti = t - index1;
		 
	  displace -= correspondences[mi][ti*4];
	  offs -= Vec3f(correspondences[mi][ti*4 + 1],
			correspondences[mi][ti*4 + 2],
			correspondences[mi][ti*4 + 3]);
	  recon += offs.dot(offs) + displace*displace;
	}
	reconRes += sqrt(recon / times.size());
	mi++;
      }
    }
  }

  reconRes /= numPoints;
  pointRes /= numPoints;

  basis->save("/tmp/ground_truth.b3d", basisRep);

  basisRep.getChannel(1).getNormalizedImage().save("/tmp/c1.png");
  basisRep.getChannel(2).getNormalizedImage().save("/tmp/c2.png");
  basisRep.getChannel(3).getNormalizedImage().save("/tmp/c3.png");

  std::cout << "Residual:" << pointRes << " recon: " << reconRes << " " << std::endl;

  return 0;
}
