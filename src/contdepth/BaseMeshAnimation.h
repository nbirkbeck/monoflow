#ifndef BASE_MESH_ANIMATION_H
#define BASE_MESH_ANIMATION_H

#include "ik/mesh.h"
#include "ik/armature.h"

#include "DisplaceUV.h"


/** \brief This is a simple class to abstract the animation of the base model.
 */
class BaseMeshAnimation {
 public:
  BaseMeshAnimation():m_p_mesh(0) { }
 
  BaseMeshAnimation(Mesh& _mesh):
   m_p_mesh(&_mesh) {}

  void setBaseGeometry(const BaseGeometry::ptr& baseGeom) {
    m_p_baseGeom = baseGeom;

    if (m_p_baseGeom) 
      m_p_baseGeom->baseMeshChanged();
  }

  void setPoseKeys(PoseKeys& keys) {
    m_poseKeys = keys;
  }

  bool loadArmature(const std::string& armatureFile) {
    if (!m_armature.read(armatureFile.c_str())) return false;
    m_armature.useEuclideanRoot();
    return true;
  }

  bool setTime(int t) {
    // lets allow animations when there is no mesh.
    if (!m_p_mesh) return true;

    if (!m_poseKeys.setCurrent(t)) 
      return false;
    m_armature.set(m_poseKeys);
    m_armature.animate(*m_p_mesh);

    if (m_p_baseGeom)
      m_p_baseGeom->baseMeshChanged();
    return true;
  }

  Armature& getArmature() {
    return m_armature;
  }

  BaseGeometry::ptr getBaseGeometry() {
    return m_p_baseGeom;
  }

 protected:
  Mesh* m_p_mesh;
  BaseGeometry::ptr m_p_baseGeom;
  Armature m_armature;
  PoseKeys m_poseKeys;
};


#endif
