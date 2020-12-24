/**
   Fit an animation to an already displaced mesh.
*/
#include <GL/gl.h>
#include "utigl/glwindow.h"
#include "ik/mesh.h"
#include "ik/armature.h"
#include "ik/MeshMotionConstraint.h"
#include <nmisc/commandline.h>
#include "BaseMeshAnimation.h"
#include "FlowBasis3D.h"
#include "utigl/fbo.h"


int main(int ac, char * av[]) {
  nacb::CommandLine cline;
  std::string geomFileName;
  std::string boneFileName;
  std::string basisFileName;
  std::string animOutFileName = "/tmp/fit.anim";
  int numFrames = 0;

  cline.registerOption("geom", "The geometry.", &geomFileName, 0);
  cline.registerOption("bones", "The input bones.", &boneFileName, 0);
  cline.registerOption("basisFile", "The basis file.", &basisFileName, 0);
  cline.registerOption("numFrames", "The number of frames.", &numFrames, 0);
  cline.registerOption("animOut", "The output file name.", &animOutFileName, 0);
  cline.parse(ac, av);

  GLWindow window(32, 32);
  FrameBufferObject fbo(1024, 1024);
  fbo.bind(1);

  Mesh mesh;
  if (!mesh.loadobj(geomFileName.c_str())) {
    std::cerr << "Cannot load mesh:" << geomFileName << std::endl;
    return -1;
  }

  Armature arm;
  if (!arm.read(boneFileName.c_str())) {
    std::cerr << "Cannot load bones:" << boneFileName << std::endl;
    return -1;
  }

  nacb::Imagef basisData;
  FlowBasis3D::ptr basis = FlowBasis3D::create(basisFileName, basisData);

  // Load in the animation (and bones)
  BaseMeshAnimation baseMeshAnimation;
  DisplacedMesh dispMesh;

  BaseGeometry::ptr baseGeometrySubdiv = 
    createBaseGeometryFromMesh<ImageBaseSkinnedGeometry>(mesh, 128, 128, dispMesh);

  PoseKeys keys(arm.getBoneNames());
  basisData = basisData.resize(128, 128);

  DisplacedMesh restMesh;

  Mesh::weight_vector_t weights;
  weights.push_back(Mesh::bone_weight_t(0, 1.0));
  dispMesh.bone_weights = Mesh::bone_weights_t(dispMesh.restVert.size(), weights);  
  arm.useEuclideanRoot();

  for (int t=0; t<numFrames; t++) {
    nacb::Imagef dispMap = basis->getDisplacementMap(basisData, t);
    dispMesh.displaceAndOffset(baseGeometrySubdiv, dispMap);

    
    if (t == 0) {
      restMesh = dispMesh;
    }


    printf("Solving %d  %ld %ld\n", t, restMesh.vert.size(), dispMesh.vert.size());
    MeshMotionConstraintSolver solver(arm, restMesh, dispMesh.vert);
    keys.append(&arm);
    printf("Solved\n");

    // This is the object.
    char fname[1024];
    snprintf(fname, 1024, "/tmp/displaced-%04d.obj", t);
    dispMesh.saveobj(fname);
  }

  keys.save(animOutFileName.c_str());
  

  return 0;
}
