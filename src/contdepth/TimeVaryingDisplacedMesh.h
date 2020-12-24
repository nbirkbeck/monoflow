#ifndef TIME_VARYING_DISPLACED_MESH
#define TIME_VARYING_DISPLACED_MESH

#include "FlowBasis3D.h"
#include "BaseMeshAnimation.h"
#include "DisplaceUV.h"


class TimeVaryingDisplacedMesh {
public:
  TimeVaryingDisplacedMesh() {
    drawBaseMesh = false;
    currentTime = -1;
    drawFlow = true;
    drawTangents = true;
    drawBones = false;
    maxTextureHeight = maxTextureWidth = 0;
  }

  void setMaximumTextureResolution(int tw, int th) {
    maxTextureWidth = tw;
    maxTextureHeight = th;
  }

  void resizeBasisData(int w, int h) {
    basisData = basisData.resize(w, h);
    baseGeometry = createBaseGeometryFromMesh<ImageBaseSkinnedGeometry>(mesh, basisData.w, basisData.h, dispMesh);
    baseMeshAnimation.setBaseGeometry(baseGeometry);
  }

  bool loadMesh(const std::string& file) {
    if (!mesh.loadobj(file.c_str())) return false;
    baseMeshAnimation = BaseMeshAnimation(mesh);
    baseMeshAnimation.setBaseGeometry(baseGeometry);
    return true;
  }

  bool setBasis(const FlowBasis3D::ptr& flowBasisIn, bool initData = true) {
    basis = flowBasisIn;
    if (initData) {
      basisData = nacb::Imagef(512, 512, basis->getNumDimensions());
      basisData = 0;
    }

    if (maxTextureWidth > 0 && maxTextureHeight > 0)
      basisData = basisData.resize(maxTextureWidth, maxTextureHeight);

    baseGeometry = createBaseGeometryFromMesh<ImageBaseSkinnedGeometry>(mesh, basisData.w, basisData.h, dispMesh);
    baseMeshAnimation.setBaseGeometry(baseGeometry);
    return true;
  }

  bool loadBasis(const std::string& basisFile) {
    basis = FlowBasis3D::create(basisFile.c_str(), basisData);
    if (!basis) {
      std::cerr << "Cannot load basis from \"" << basisFile << "\"" << std::endl;
      return false;
    }
    // Set the basis, but don't init the basisData
    return setBasis(basis, false);
  }

  void setTimePoseOnly(double t) {
     baseMeshAnimation.setTime((int)round(t));
  }

  bool setTime(double t) {
    if ((int)round(t) != currentTime) {
      currentTime = round(t);
      baseMeshAnimation.setTime((int)round(t));

      if (basis) {
	nacb::Imagef dispMap = basis->getDisplacementMap(basisData, (int)round(t));
	dispMesh.displaceAndOffset(baseGeometry, dispMap);
      }
    }
    return true;
  }

  bool loadKinematics(const std::string& bonesFileName,
		      const std::string& animFileName) 
  {
    PoseKeys poseKeys;
    baseMeshAnimation = BaseMeshAnimation(mesh);
      
    if (!poseKeys.load(animFileName.c_str())) {
      std::cerr << "Error loading animation from:" << animFileName << std::endl;
      return false;
    }
    baseMeshAnimation.setPoseKeys(poseKeys);

    if (!baseMeshAnimation.loadArmature(bonesFileName.c_str())) {
      std::cerr << "Error loading armature from " << bonesFileName << std::endl;
      return false;
    }

    baseMeshAnimation.setBaseGeometry(baseGeometry);
    return true;
  }

  void draw() { 
    if (drawBaseMesh) {
      glPushAttrib(GL_ALL_ATTRIB_BITS);
      glLineWidth(2);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      mesh.draw();

      glPopAttrib();
    }
    else
      dispMesh.draw();

    if (baseMeshAnimation.getArmature().root && drawBones)
      baseMeshAnimation.getArmature().root->draw();

    if (drawFlow && baseGeometry) {
      glDisable(GL_LIGHTING);
      glDisable(GL_TEXTURE_2D);

      glColor3f(0.6, 0.0, 0.0);
      std::vector<std::vector<nacb::Vec3f> > pointsByTime((basisData.w/4 + 1) * (basisData.h/4 + 1));

      // Compute the animated points so we can draw lines.
      for (double t = max(0.0, currentTime - 2); t <= currentTime + 10; t++) {
	baseMeshAnimation.setTime(int(round(t)));
	
	int di = 0;
	for (int y=0; y<basisData.h; y+=4) {
	  for (int x=0; x<basisData.w; x+=4, di++) {
	    double displace = 0.0;
	    nacb::Vec3f offset = (*basis)(displace, x, y, basisData, int(round(t)));
	    nacb::Vec3f point = baseGeometry->displaceAndOffset(x, y, displace, offset.x, offset.y, offset.z);
	    pointsByTime[di].push_back(point);
	  }
	}
      }
      glDisable(GL_LIGHTING);
      for (int i=0; i<(int)pointsByTime.size(); i++) {
	glBegin(GL_LINE_STRIP);
	for (int j=0; j<(int)pointsByTime[i].size(); j++) {
	  glVertex3fv(pointsByTime[i][j].data);
	}
	glEnd();
      }
    }
    if (drawTangents) {
      glBegin(GL_LINES);
      for (int y=0; y<basisData.h; y+=4) {
	for (int x=0; x<basisData.w; x+=4) {
	  nacb::Vec3f tangents[3], normal;
	  nacb::Vec3f co = baseGeometry->displaceAndOffset(x, y, 0.0, 0, 0, 0, &normal, tangents);
	  
	  glColor3f(1, 0, 0);
	  nacb::Vec3f v = co + tangents[0]*0.05;
	  glVertex3fv(co.data);
	  glVertex3fv(v.data);

	  glColor3f(0, 1, 0);
	  v = co + tangents[1]*0.05;
	  glVertex3fv(co.data);
	  glVertex3fv(v.data);

	  glColor3f(0, 0, 1);
	  v = co + tangents[2]*0.05;
	  glVertex3fv(co.data);
	  glVertex3fv(v.data);
	}
      }
      glEnd();
    }
    baseMeshAnimation.setTime(int(round(currentTime)));
  }

  BaseGeometry::ptr getBaseGeometry() {
    return baseGeometry;
  }

public:
  bool drawBaseMesh;
  bool drawFlow;
  bool drawTangents;
  bool drawBones;
  double currentTime;
  Mesh mesh;
  BaseMeshAnimation baseMeshAnimation;	       
  FlowBasis3D::ptr basis;
  nacb::Imagef basisData;

  BaseGeometry::ptr baseGeometry;
  DisplacedMesh dispMesh;

  int maxTextureWidth, maxTextureHeight;
};


#endif 
