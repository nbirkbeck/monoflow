//#pragma implementation "armature.h"
#include "armature.h"

void PoseKeys::append(Armature * arm){
  std::vector<Vec3d>       pos;
  std::vector<Quaternion> quat;
  std::vector<int> withPos;
  
  for(int i=0; i<(int)arm->size(); i++){
    Bone * bone = (*arm)[i];

    // Always add the pose for the root.
    if((bone->joint && bone->joint->useloc) || !bone->parent) {
      pos.push_back(bone->getPoseLoc());
      withPos.push_back(i);
    }
    
    quat.push_back(bone->getPoseQuat());
  }
  append(quat, pos, withPos);
}

void Armature::relaxHips(bool test){
  Bone * back=operator[]("lowerback");
  Bone * lfemur=operator[]("lfemur");
  Bone * rfemur=operator[]("rfemur");
  Bone * thorax=operator[]("thorax");

  if(!back || !lfemur || !rfemur || !thorax)return;

  //Hack for upper back parts
  Bone * finger=thorax;
  while(finger->name!="lowerback"){
    if(3==finger->getNParams()){
      Bone * parent = finger->parent;
      while(parent && 3!=parent->getNParams())
	parent=parent->parent;
      if(parent){
	Mat4x4 cur=finger->tform;
	
	double parms[3];
	finger->getParams(parms);
	parms[1]*=0.5;
	finger->setParams(parms,false,true);

	Mat4x4 fixed=(parent->tform*parent->lanim.eucInv())*parent->lmat;
	Mat4x4 panim=fixed.eucInv()*cur*finger->lanim.eucInv()*Mat4x4::trans(0,-parent->len,0);
	parent->setPoseQuat(Quaternion::fromMatrix(panim));
	cout << "Panim:" << panim << endl;
      }
    }
    finger=finger->parent;
  }
  
  

  double parms[3];
  back->getParams(parms);
  Mat4x4 cur=back->tform;//.copy();
      
  Mat4x4 lcur=lfemur->tform;
  Mat4x4 rcur=rfemur->tform;

  if(test && fabs(parms[1])<1e-4){
    parms[1]=(double(rand())/double(RAND_MAX))*2.0-1.0;
    printf("setting y to %f\n",parms[1]);
  }
  else
    parms[1]*=0.5;
  back->setParams(parms,false,true);

  Mat4x4 del=cur*back->tform.eucInv();
  //lmat*(lmat.inv()*upd*lmat)*parms*trans*(child)=

  Mat4x4 rnew=root->lmat.eucInv()*del*root->lmat;
  root->joint->apply(rnew);
  root->setPoseLoc(Vec3d(rnew(0,3),rnew(1,3),rnew(2,3)));
  root->setPoseQuat(Quaternion::fromMatrix(rnew));
  root->updateTransform(true);

  Mat4x4 laft=lfemur->tform;
  Mat4x4 raft=rfemur->tform;
      
  Mat4x4 tmp=((raft*rfemur->lanim.eucInv())*rfemur->lmat).eucInv()*rcur;
  cout << "q:" << Quaternion::fromMatrix(tmp) << endl;
  cout << "was:" << rfemur->getPoseQuat() << endl;

  Quaternion qres=Quaternion::slerp(0.5,Quaternion::fromMatrix(tmp),
				    rfemur->getPoseQuat());
  qres.normalize();
  rfemur->setPoseQuat(qres);
  
  tmp=((laft*lfemur->lanim.eucInv())*lfemur->lmat).eucInv()*lcur;
  qres=Quaternion::slerp(0.5,Quaternion::fromMatrix(tmp),
			 lfemur->getPoseQuat());
  qres.normalize();
  lfemur->setPoseQuat(qres);

  //Just reduce the y-rotation on femurs
  lfemur->getParams(parms);
  parms[1]*=0.5;
  lfemur->setParams(parms,false,true);

  rfemur->getParams(parms);
  parms[1]*=0.5;
  rfemur->setParams(parms,false,true);
  //Done with the femurs

  synchMesh();
  //animate(globs.mesh);
}


double * Armature::meshJacobian(Mesh & mesh,
				const double * animp,
				int animp_len,
				bool differential,
				JacobianCache * jout){
  double * back=new double[animp_len];
  root->getParamsAll(back);

  int nvert=mesh.vert.size();

  //The tforms have to be current but the mesh.vert doesn't
  root->setParamsAll(animp,differential);
  root->updateTransform(true);

  JacobianCache * jcache=jout;
  if(!jcache)jcache=new JacobianCache();
  if(jcache->nvert!=nvert ||
     jcache->nbone!=(int)boneByID.size() ||
     jcache->nparm!=animp_len){

    jcache->alloc(nvert,animp_len,boneByID.size());

    //Each vertex is: \sum_j (B_j^t(\alpha) (B_j^0)^-1 x_i) w_j
    //Need to find the parameters that affect each bone transformation first:
    //   The bone tform is affected by it's own parameters and that of all anscestors
    vector<int> * bone_affected_by=new vector<int>[boneByID.size()];
    
    root->getParamsFromAll(jcache->from,jcache->index);

    //make list of parameters that came from each bone
    vector<int> * bone_param_inds=new vector<int>[boneByID.size()];
    for(int i=0;i<animp_len;i++)bone_param_inds[jcache->from[i]].push_back(i);
    
    //Each bone is affected by the parameters that affect its parent
    //and its own parameters.  This is okay as indices are breadth first
    for(int i=0;i<(int)boneByID.size();i++){
      Bone * bone=boneByID[i];
      if(bone->parent)
	bone_affected_by[i].insert(bone_affected_by[i].end(),
				   bone_affected_by[bone->parent->id].begin(),
				   bone_affected_by[bone->parent->id].end());
      bone_affected_by[i].insert(bone_affected_by[i].end(),
				 bone_param_inds[i].begin(),
				 bone_param_inds[i].end());
      printf("%s is affected by %ld\n",
             bone->name.c_str(),bone_affected_by[i].size());

      //Compute the parameters that affect bone i 
      for(int k=0;k<(int)bone_affected_by[i].size();k++){
	jcache->bone_affected[i*animp_len+bone_affected_by[i][k]]=true;
      }
    }
    bool * vert_affected_by=new bool[animp_len];
    bool * vert_affected_by_bone=new bool[boneByID.size()];

    //Compute the vertices affected by bone i 
    for(int i=0;i<(int)mesh.restVert.size();i++){
      memset(vert_affected_by,0,sizeof(bool)*animp_len);
      
      for(int k=0;k<(int)mesh.bone_weights[i].size();k++){
	int bi=mesh.bone_weights[i][k].bone;
	for(int h=0;h<animp_len;h++)
	  vert_affected_by[h]|=jcache->bone_affected[bi*animp_len+h];
      }
      int cnt=0;
      for(int h=0;h<animp_len;h++)
	cnt+=vert_affected_by[h];
      
      //A vertex is affected by a bone if any of the parameters that affect
      //the bone also affect the vertex.
      memset(vert_affected_by_bone,0,sizeof(bool)*boneByID.size());
      for(int h=0;h<animp_len;h++)
	if(vert_affected_by[h])
	  vert_affected_by_bone[jcache->from[h]]=true;

      int cnt2=0;
      for(int h=0;h<(int)boneByID.size();h++)
	if(vert_affected_by_bone[h]){
	  cnt2++;
	  jcache->affected_vert[h].push_back(i);
	}
    }
    delete [] vert_affected_by;
    delete [] vert_affected_by_bone;

    for(int h=boneByID.size(); h < (int)boneByID.size(); h++){
      printf("bone %s affects %ld/%ld verts\n",
             boneByID[h]->name.c_str(),jcache->affected_vert[h].size(),mesh.restVert.size());
    }
    delete [] bone_param_inds;
    delete [] bone_affected_by;
  }
  double * deriv=new double[nvert*3*animp_len];

  StopWatch watch(true);
  for(int pi=0;pi<animp_len;pi++){
    Bone * bone=boneByID[jcache->from[pi]];
    bone->setDeriv(jcache->index[pi],differential);

    for(int k=0;k<(int)bone->kids.size();k++)
      bone->kids[k]->updateTransform(true);
	
    double * col=deriv+pi*3*nvert;
    memset(col,0,nvert*3*sizeof(double));

    for(int i=0;i<(int)jcache->affected_vert[bone->id].size();i++){
      int vi=jcache->affected_vert[bone->id][i];
      Vec4d v(mesh.restVert[vi].x,
	      mesh.restVert[vi].y,
	      mesh.restVert[vi].z,1);
      for(int k=0; k<(int)mesh.bone_weights[vi].size(); k++){
	Mesh::bone_weight_t & bw=mesh.bone_weights[vi][k];

	//Is this vertex affected by the current parameter
	if(jcache->bone_affected[bw.bone*animp_len+pi])
	{
	  Vec4d tp=boneByID[bw.bone]->tform*(boneByID[bw.bone]->irest*v);
	  
	  col[3*vi  ]+=tp.x*bw.weight;
	  col[3*vi+1]+=tp.y*bw.weight;
	  col[3*vi+2]+=tp.z*bw.weight;
	}
      }
    }
    //Reset the transform to be the current value
    bone->updateTransform(true);
  }
  
  //printf("loop:%f \n",double(watch));
  root->setParamsAll(back);
  root->updateTransform(true);

  delete [] back;
  //Using self made cache, delete it
  if(!jout){
    fprintf(stderr,"Deleting jacobian cache\n");
    delete jcache;
  }
  return deriv;
}


double * Armature::meshJacobian_fd(Mesh & mesh,
				   const double * animp_const,
				   int animp_len,bool differential){
  double * animp=new double[animp_len];
  memcpy(animp,animp_const,sizeof(double)*animp_len);

  double * back=new double[animp_len];
  root->getParamsAll(back);

  root->setParamsAll(animp,differential);
  root->updateTransform(true);

  Mat4x4 * tforms=new Mat4x4[boneByID.size()];
  for(int i=0;i<(int)boneByID.size();i++){
    tforms[i]=boneByID[i]->tform*boneByID[i]->irest;
  }
  mesh.animate(tforms);
  vector<Vec3f> vcur=mesh.vert;

  //Test the derivative of the vertices w.r.t to the parameters
  double * deriv_fd=new double[vcur.size()*3*animp_len];

  for(int pi=0;pi<animp_len;pi++){
    double c=animp[pi];
    double dx=1e-2;
    animp[pi]=c+dx;
    root->setParamsAll(back);
    root->setParamsAll(animp,differential);
    root->updateTransform(true);
    for(int i=0;i<(int)boneByID.size();i++){
      tforms[i]=boneByID[i]->tform*boneByID[i]->irest;
    }
    mesh.animate(tforms);
    vector<Vec3f> vup=mesh.vert;

    animp[pi]=c-dx;
    root->setParamsAll(back);
    root->setParamsAll(animp,differential);
    root->updateTransform(true);
    for(int i=0;i<(int)boneByID.size();i++){
      tforms[i]=boneByID[i]->tform*boneByID[i]->irest;
    }
    mesh.animate(tforms);
    vector<Vec3f> vdn=mesh.vert;
    
    animp[pi]=c;
    
    double * col=deriv_fd+pi*3*vcur.size();
    for(int i=0;i<(int)vcur.size();i++){
      Vec3f d=(vup[i]-vdn[i])*(1.0/(2.0*dx));
      col[3*i  ]=d.x;
      col[3*i+1]=d.y;
      col[3*i+2]=d.z;
    }
  }
  mesh.vert=vcur;//restore the current verts

  root->setParamsAll(back);
  root->updateTransform(true);

  delete [] animp;
  
  delete [] tforms;
  delete [] back;
  return deriv_fd;
}


void Armature::applyCollisionConstraint(Mesh & mesh){
  if (!mesh.spheres.size()) return;

  vector<pair<int,int> > cons;
  vector<double> dists;

  int nparms=root->getNParamsAll();
  int * from = new int[nparms];
  int * index = new int[nparms];

  //Each vertex is: \sum_j (B_j^t(\alpha) (B_j^0)^-1 x_i) w_j
  //Need to find the parameters that affect each bone transformation first:
  //   The bone tform is affected by it's own parameters and that of all anscestors
  vector<int> * bone_affected_by=new vector<int>[boneByID.size()];
    
  root->getParamsFromAll(from,index);

  //make list of parameters that came from each bone
  vector<int> * bone_param_inds=new vector<int>[boneByID.size()];
  for(int i=0;i<(int)nparms;i++)bone_param_inds[from[i]].push_back(i);

  bool * bone_affected = new bool[boneByID.size()*nparms];
  memset(bone_affected, 0, sizeof(bool)*boneByID.size()*nparms);

  //Each bone is affected by the parameters that affect its parent
  //and its own parameters.  This is okay as indices are breadth first
  for(int i=0;i<(int)boneByID.size();i++){
    Bone * bone=boneByID[i];
    if(bone->parent)
      bone_affected_by[i].insert(bone_affected_by[i].end(),
				 bone_affected_by[bone->parent->id].begin(),
				 bone_affected_by[bone->parent->id].end());
    bone_affected_by[i].insert(bone_affected_by[i].end(),
			       bone_param_inds[i].begin(),
			       bone_param_inds[i].end());
    //printf("collision: %s is affected by %d\n",bone->name.c_str(),bone_affected_by[i].size());
    
    //Compute the parameters that affect bone i 
    for(int k=0;k<(int)bone_affected_by[i].size();k++){
      bone_affected[i*nparms+bone_affected_by[i][k]]=true;
    }
  }
  vector<int> * bone_affects = new vector<int>[boneByID.size()];

  int nspheres=mesh.spheres.size();
  bool * sphere_affected_by = new bool[nparms];

  for(int i=0;i<nspheres;i++){
    memset(sphere_affected_by, 0, nparms);
    for(int j=0;j<(int)mesh.spheres[i].weights.size();j++){
      int bi=mesh.spheres[i].weights[j].bone;
      for(int k=0;k<(int)bone_affected_by[bi].size();k++){
	sphere_affected_by[bone_affected_by[bi][k]]=true;
      }
    }
    for(int k=0;k<nparms;k++){
      if(sphere_affected_by[k]){
	int bi=from[k];
	if(find(bone_affects[bi].begin(),bone_affects[bi].end(),i)
	   ==bone_affects[bi].end())
	  bone_affects[bi].push_back(i);
      }
    }
  }
  for(int i=0;i<(int)boneByID.size();i++){
    std::sort(bone_affects[i].begin(),bone_affects[i].end());
  }

  delete [] sphere_affected_by;

  double * pback = new double[nparms];
  double * pcur = new double[nparms];
  int * sphere_map = new int[nspheres];
  
  root->getParamsAll(pback);
  memcpy(pcur,pback,sizeof(double)*nparms);

  for(int its=0;its<20;its++){
    cons.clear();
    dists.clear();

    int used_spheres=0;
    memset(sphere_map,0xFF,sizeof(int)*nspheres);

    for(int i=0;i<(int)mesh.spheres.size();i++){
      for(int j=0;j<(int)mesh.spheres[i].interacts.size();j++){
	int other=mesh.spheres[i].interacts[j];
	if(i<other){
	  Vec3f diff=(mesh.spheres[i].c-mesh.spheres[other].c);
	  double d=diff.len();
	  double tr=(mesh.spheres[i].r+mesh.spheres[other].r)*0.9;
	  if(d<=tr){
	    //printf("got collision between %d and %d  %f\n",i,other,tr-d);
	    //(theta_i-curtheta_i)^2
	    //(mesh.spheres[i].c-mesh.spheres[other].c)^2-r^2
	    cons.push_back(pair<int,int>(i,other));
	    dists.push_back((tr+1e-4)*(1.01));

	    if(sphere_map[i]<0)sphere_map[i]=used_spheres++;
	    if(sphere_map[other]<0)sphere_map[other]=used_spheres++;
	  }
	}
      }
    }
    Matrix J(nparms+cons.size(),nparms);
    Matrix f(nparms+cons.size(),1);
    J.setAll(0);
    for(int i=0;i<nparms;i++){
      J(i,i)=0.1;
      f[i]=-(pcur[i]-pback[i])*0.1;
    }
    vector<int> naffects;
    for(int i=0;i<(int)boneByID.size();i++){
      int withcons=0;
      for(int j=0;j<(int)bone_affects[i].size();j++){
	withcons+=(sphere_map[bone_affects[i][j]]>=0);
      }
      naffects.push_back(withcons);
    }
    Matrix SJ(3*used_spheres,nparms);
    SJ.setAll(0);

    for(int pi=0;pi<nparms;pi++){
      int bi=from[pi];

      if(naffects[bi]){
	Bone * bone = boneByID[bi];
	bone->setDeriv(index[pi],false);
      
	for(int k=0;k<(int)bone->kids.size();k++)
	  bone->kids[k]->updateTransform(true);

	for(int j=0;j<(int)bone_affects[bi].size();j++){
	  int si=bone_affects[bi][j];
	  if(sphere_map[bone_affects[bi][j]]>=0){
	    Vec4d v(mesh.spheres[si].restc.x,
		    mesh.spheres[si].restc.y,
		    mesh.spheres[si].restc.z,1);
	  
	    for(int k=0;k<(int)mesh.spheres[si].weights.size();k++){
	      Mesh::bone_weight_t & bw=mesh.spheres[si].weights[k];

	      if(bone_affected[bw.bone*nparms+pi]){
		Vec4d tp=boneByID[bw.bone]->tform*(boneByID[bw.bone]->irest*v);
	      
		int ui=sphere_map[si];
		SJ(3*ui  ,pi)+=tp.x*bw.weight;
		SJ(3*ui+1,pi)+=tp.y*bw.weight;
		SJ(3*ui+2,pi)+=tp.z*bw.weight;
	      }
	    }
	  }
	}
	//Reset the transform to be the current value
	bone->updateTransform(true);
      }
    }

    for(int i=0;i<(int)cons.size();i++){
      const int s1=cons[i].first;
      const int s2=cons[i].second;
      int m1=sphere_map[s1];
      int m2=sphere_map[s2];
      //norm(p1-p2)=r+0.01
      Vec3f d=(mesh.spheres[s1].c-mesh.spheres[s2].c);
      double len=d.len();
      for(int j=0;j<nparms;j++){
	double x=SJ(3*m1  ,j)-SJ(3*m2,  j);
	double y=SJ(3*m1+1,j)-SJ(3*m2+1,j);
	double z=SJ(3*m1+2,j)-SJ(3*m2+2,j);

	J(nparms+i,j)=(d.x*x+d.y*y+d.z*z)/(len);
      }
      f[nparms+i]=-std::min(len-dists[i],0.0);
    }
    Matrix upd=Matrix::LlinLeastSq(J,f);
    //upd.transpose().printMatlab("upd");
    printf("residual before : %f\n",f.dot(f));
  
    for(int i=0;i<nparms;i++){
      pcur[i]+=upd[i];
    }
    root->setParamsAll(pcur,false);
    animate(mesh);
    root->getParamsAll(pcur,false);
  }

  delete [] sphere_map;
  delete [] pback;

  delete [] bone_affects;
  delete [] bone_affected;

  delete [] bone_param_inds;
  delete [] bone_affected_by;
  delete [] from;
  delete [] index;
}


void Armature::initHumanTips(){
  vector<string> tip_names;
  vector<int> lens;
  
  tip_names.push_back("lradius");
  tip_names.push_back("rradius");
  
  tip_names.push_back("ltibia");
  tip_names.push_back("rtibia");
  tip_names.push_back("lhumerus");
  tip_names.push_back("rhumerus");
  
  tip_names.push_back("lfemur");
  tip_names.push_back("rfemur");
  
  
  lens.push_back(2);
  lens.push_back(2);
  lens.push_back(2);
  lens.push_back(2);
  lens.push_back(1);
  lens.push_back(1);
  
  lens.push_back(1);
  lens.push_back(1);
  

  tip_names.push_back("left_shoulder");
  tip_names.push_back("left_elbow");
  tip_names.push_back("left_wrist");
  tip_names.push_back("left_hand");
  lens.push_back(1);
  lens.push_back(2);
  lens.push_back(3);
  lens.push_back(4);


  tip_names.push_back("right_shoulder");
  tip_names.push_back("right_elbow");
  tip_names.push_back("right_wrist");
  tip_names.push_back("right_hand");
  lens.push_back(1);
  lens.push_back(2);
  lens.push_back(3);
  lens.push_back(4);

  tip_names.push_back("left_hip");
  tip_names.push_back("left_knee");
  tip_names.push_back("left_ankle");
  tip_names.push_back("left_toe");
  lens.push_back(1);
  lens.push_back(2);
  lens.push_back(3);
  lens.push_back(4);

  tip_names.push_back("right_hip");
  tip_names.push_back("right_knee");
  tip_names.push_back("right_ankle");
  tip_names.push_back("right_toe");
  lens.push_back(1);
  lens.push_back(2);
  lens.push_back(3);
  lens.push_back(4);

  tip_names.push_back("neck");
  tip_names.push_back("head");
  tip_names.push_back("head_right");
  tip_names.push_back("head_left");
  tip_names.push_back("head_back");
  tip_names.push_back("head_up");

  lens.push_back(1);
  lens.push_back(2);
  lens.push_back(3);
  lens.push_back(3);
  lens.push_back(3);
  lens.push_back(3);

  setIKDepths(tip_names, lens, true);
}


std::pair<double, Matrix>  compareParamsDirect(Armature & arm, const Matrix & p1, const Matrix & p2){
  Matrix diff;
  diff = p1 - p2;

  double avg = 0.0;
  for(int i=0; i<diff.m*diff.n; i++){
    diff[i] *= diff[i];
    avg += diff[i];
  }

  return std::pair<double, Matrix>(avg/(diff.m*diff.n), diff);
}


std::pair<double, Matrix> compareParamsQuats(Armature & arm, const Matrix & p1, const Matrix & p2){
  //FIXME: euclidean other bones.
  Matrix res(arm.size()+1, 1);
  
  arm.setParams(p1.data);
  arm.updateTransform();

  Vec3d l1 = arm.root->getPoseLoc();
  std::vector<Quaternion> quats;
  for(int i=0; i<(int)arm.size(); i++){
    quats.push_back(arm[i]->getPoseQuat());
  }

  arm.setParams(p2.data);
  arm.updateTransform();

  double avg = 0;
  res[0] = (l1 - arm.root->getPoseLoc()).len();
  res[0] = res[0]*res[0];
  
  avg += res[0];

  for(int i=0; i<(int)arm.size(); i++){
    Quaternion q2 = arm[i]->getPoseQuat();
    double d = q2.dot(quats[i]);
    if(d<0)d = -d;
    
    res[i+1] = acos(std::min(1.0, d));
    avg += res[i+1];
  }

  return std::pair<double, Matrix>(avg/(arm.size()+1), res);
}


void compareParamsSummary(Armature & arm, PoseKeys & k1, PoseKeys & k2){
  int nframes = std::min(k1.getNumFrames(), k2.getNumFrames());
  int nparms = arm.getNParams();
  Matrix p1(nparms, 1), p2(nparms, 1);

  Matrix d(nframes, 1), q(nframes, 1);

  for(int i=0; i<nframes; i++){
    k1.setCurrent(i);
    arm.set(k1);
    arm.updateTransform();
    arm.getParams(p1.data);

    k2.setCurrent(i);
    arm.set(k2);
    arm.updateTransform();
    arm.getParams(p2.data);

    std::pair<double, Matrix> d_err = compareParamsDirect(arm, p1, p2);
    std::pair<double, Matrix> q_err = compareParamsQuats(arm, p1, p2);

    d[i] = d_err.first;
    q[i] = q_err.first;

    printf("d_err(%d) = %lf\n", i, d_err.first);
    printf("q_err(%d) = %lf\n", i, q_err.first);
  }
  //This was not so meaningful for small errors.
  //printf("Summary: %lf %lf\n", d.dot(d), q.dot(q));
  printf("Summary: %lf %lf\n", sqrt(d.dot(d)/nframes), sqrt(q.dot(q)/nframes));
}


