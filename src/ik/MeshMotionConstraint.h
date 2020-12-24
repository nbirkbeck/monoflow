#ifndef MESH_MOTION_CONSTRAINT_H
#define MESH_MOTION_CONSTRAINT_H

#include "armature.h"
#include "mesh.h"
#include <nmath/vec3.h>
#include <vector>


class MeshMotionConstraintSolver {
public:
  MeshMotionConstraintSolver(Armature& _arm,
			     Mesh & _mesh,
			     nacb::Vec3d * _pcons,
			     int * _pcnt): arm(_arm), mesh(_mesh), pcons(_pcons), pcnt(_pcnt), printRes(true) {
    for (int i=0; i<(int)mesh.vert.size(); i++) {
      if (pcnt[i]) {
	withCons.push_back(i);
      }
    }
    solve();
  }

  MeshMotionConstraintSolver(Armature& _arm,
			     Mesh & _mesh,
			     const std::vector<nacb::Vec3f>& _cons): arm(_arm), mesh(_mesh), pcons(0), pcnt(0), printRes(false) {
    std::vector<nacb::Vec3d> pconsd(_cons.size());
    pcons = &(pconsd[0]);

    for (int i=0; i<(int)mesh.vert.size(); i++) {
      pconsd[i] = nacb::Vec3d(_cons[i].x,
			      _cons[i].y,
			      _cons[i].z);			      
      withCons.push_back(i);
    }
    solve();
  }

  void solve() {
    nacb::Matrix fvec(withCons.size()*3, 1);
    nacb::Matrix x(arm.getNParams(), 1);

    arm.getParams(x.data);

    // Only do a few iterations here (no need optimizing completely
    nacb::Matrix::lsqnonlin(s_objectiveFunc, this, x, fvec, 1e-6, 1e-6, 1e-6, 12);

    arm.setParams(x.data);
    arm.animate(mesh);
  }

  static int s_objectiveFunc(void * p,
			     int m,
			     int n,
			     const double * x,
			     double * fvec,
			     double * fjac,
			     int ldfjac,
			     int iflag)
  {
    return ((MeshMotionConstraintSolver*)p)->objectiveFunc(p, m, n, x, fvec, fjac, ldfjac, iflag);
  }

  int objectiveFunc(void * p,
		    int m,
		    int n,
		    const double * x,
		    double * fvec,
		    double * J,
		    int ldfjac,
		    int iflag) {
    arm.setParams(x);
    arm.animate(mesh);
    
    if (iflag == 1) {
      double res = 0;
      for (int i=0; i<(int)withCons.size(); i++) {
	int vi = withCons[i];

	fvec[3*i + 0] = mesh.vert[vi].x - pcons[vi].x;
	fvec[3*i + 1] = mesh.vert[vi].y - pcons[vi].y;
	fvec[3*i + 2] = mesh.vert[vi].z - pcons[vi].z;

	res += fvec[3*i + 0]*fvec[3*i + 0];
	res += fvec[3*i + 1]*fvec[3*i + 1];
	res += fvec[3*i + 2]*fvec[3*i + 2];
      }
      if (printRes)
	printf("%f\n", res);
    }
    else if (iflag == 2) {
      static JacobianCache jcache;
      double * jac=arm.meshJacobian(mesh, x, n, false, &jcache);

      for(int i=0; i<(int)withCons.size(); i++){
	for(int j=0; j<n; j++){
	  J[3*i + 0 + j*ldfjac] =jac[3*withCons[i] + 0 + j*3*mesh.vert.size()];
	  J[3*i + 1 + j*ldfjac] =jac[3*withCons[i] + 1 + j*3*mesh.vert.size()];
	  J[3*i + 2 + j*ldfjac] =jac[3*withCons[i] + 2 + j*3*mesh.vert.size()];
	}
      }
      delete [] jac;
    }
    return iflag;
  }

  nacb::Vec3d * pcons;
  int * pcnt;
  Armature& arm;
  Mesh& mesh;
  std::vector<int> withCons;
  bool printRes;
};


// Old version (simple, only did gauss-newton styl eminimization)
void solveMeshMotionConstraintSimple(Armature & arm,Mesh & mesh,Vec3d * pcons,int * pcnt,int maxits);

#endif
