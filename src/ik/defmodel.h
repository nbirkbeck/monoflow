/** \class DeformationModel
    
    This class represents how a mesh deforms for a given frame, or pose
    of an animation.  The deform function takes the frame parameter (instead
    of just the transformations) so that the deformation can be based on time).

    Issues: 
    -Instead of having a deform function, should the DeformationModel
     have references to the mesh, armature and keys?

    -Raw displacements and pose keys are key-framed but do not keep track
     of which frame (in a sequence) they are from, implying they must have
     the same skip (and offset).

     Warning:  Make sure the armature transforms are up-to-date (e.g., arm.updateTransform()).
               This is done in both the arm.animate, and arm.getTransforms, so there is
	       probably no need to worry.
*/
#ifndef DEFMODEL_H
#define DEFMODEL_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "displace.h"
#include "armature.h"
#include "mesh.h"


class DeformationModel {
public:

  class ConstructorInfo {
  public:
    const Mesh     * mesh;
    const Armature  * arm;
    std::string  filename;
    
    ConstructorInfo(const Mesh * _mesh = 0,
		    const Armature * _arm = 0,
		    const std::string & _filename = ""):
    mesh(_mesh), arm(_arm), filename(_filename)
    {
       
    }

    ConstructorInfo(const Mesh & _mesh,
		    const Armature & _arm,
		    const std::string & _filename = ""):
    mesh(&_mesh), arm(&_arm), filename(_filename)
    {
       
    }
  };

  DeformationModel(){ }
  
  virtual ~DeformationModel(){ }

  virtual bool deform(Mesh & mesh, Armature & arm, int frame){
    arm.animate(mesh);
    return true;
  }

  //Convenience for single line deformation
  virtual bool deform(Mesh & mesh, Armature & arm, PoseKeys & keys){
    arm.set(keys);
    return deform(mesh, arm, keys.getCurrent());
  }

  static DeformationModel * create(const DeformationModel::ConstructorInfo & cinfo){
    return new DeformationModel();
  }
};


typedef boost::shared_ptr<DeformationModel> DeformationModelPtr;


class RawDeformationModel: public DeformationModel {
public:  
  RawDeformationModel(){ }
  
  RawDeformationModel(const Mesh & mesh, 
		      const std::string & fname){
    disps =  dpl::read_dpl(fname.c_str(), mesh);
  }

  RawDeformationModel(const std::vector<dpl::displace_t> & _disps){
    disps = _disps;
  }
  
  virtual ~RawDeformationModel(){ }

  bool deform(Mesh & mesh, Armature & arm, int frame){
    if(frame>=disps.size())return false;    

    std::vector<Mat4x4> transforms = arm.getTransforms();
    for(int i=0; i<mesh.vert.size(); i++){
      mesh.vert[i] = dpl::transform(transforms, mesh.bone_weights[i], 
				    mesh.restVert[i] + disps[frame][i], 1.0);
    }    
    return true;
  }

  static DeformationModel * create(const ConstructorInfo & cinfo){
    return new RawDeformationModel(*cinfo.mesh, cinfo.filename);
  }

protected:
  std::vector<dpl::displace_t> disps;
};


class WorldObjDeformationModel: public DeformationModel {
public:  
  WorldObjDeformationModel(){ }
  WorldObjDeformationModel(const std::string & _fname) : fname(_fname){ }

  bool deform(Mesh & mesh, Armature & arm, int frame){
    Mesh dup;
    dup.loadobj(nstr::printf(fname.c_str(), frame).c_str());

    assert(dup.vert.size() == mesh.vert.size());
    mesh.vert = dup.vert;

    std::vector<Mat4x4> transforms = arm.getTransforms();
    return true;
  }

  static DeformationModel * create(const ConstructorInfo & cinfo){
    return new WorldObjDeformationModel(cinfo.filename);
  }

protected:
  const std::string fname;
};

class DeformationModelFactory {
 public:
  typedef DeformationModel * (* creator_func_t)(const DeformationModel::ConstructorInfo &);
    
  static std::map<std::string,  creator_func_t> & getCreators(){
    static std::map<std::string, creator_func_t> creators;
    
    if(creators.size()==0){
      creators["basic"] = &DeformationModel::create;
      creators["dpl"] = &RawDeformationModel::create;
      creators["world"] = &WorldObjDeformationModel::create;
    }
    return creators;
  }

  static DeformationModelPtr create(const std::string & type,
				    const DeformationModel::ConstructorInfo & cinfo){
    if(getCreators().count(type)){
      return DeformationModelPtr(getCreators()[type](cinfo));
    }
    return DeformationModelPtr(new DeformationModel());
  }
};


#endif
