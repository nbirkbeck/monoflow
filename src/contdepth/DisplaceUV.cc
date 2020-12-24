#include "DisplaceUV.h"
#include "ik/triangle.h"
#include "utigl/glcommon.h"

using namespace std;
using namespace nacb;

/**
   Create the image-based representation of the mesh.
   
   Also create a full resolution subdivided mesh that can be displaced
*/
template <class BaseGeometryType>
BaseGeometry::ptr createBaseGeometryFromMesh(Mesh & mesh, int w, int h, 
					     DisplacedMesh & submesh){
  nacb::Imagef uv_co(w, h, 3), uv_no(w, h, 3);
  nacb::Imagef mask(w, h, 1);

  nacb::Image32 timage(w, h, 1);
  nacb::Image32 vertMap(w, h, 1);
  nacb::Imagef baryImage = getBaryImage(mesh, w, h, true, &timage);

  std::vector<Triangle> triangles;
  for (int i=0; i<(int)mesh.tris.size(); i++) {
    Vec3f v[3] = {mesh.tvert[mesh.tris[i].tci[0]],
		  mesh.tvert[mesh.tris[i].tci[1]],
		  mesh.tvert[mesh.tris[i].tci[2]]};
    v[0].z = v[1].z = v[2].z = 0;
    triangles.push_back(Triangle(v[0], v[1], v[2]));
  }

  mask = 0;
  uv_co = 0;
  uv_no = 0;

  vertMap = 0;
  submesh = DisplacedMesh();

  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      int trii = timage(x, y);

      if (trii < (int)mesh.tris.size()) {
	Vec3f tc(double(x)/w, double(y)/h, 0.0);
	Vec3f dir(0.f, 0.f, 1.f);
	
	Vec3f bary(baryImage(x, y, 0),
		   baryImage(x, y, 1),
		   baryImage(x, y, 2));
		
	double t, a, b, c;
	if (triangles[trii].intersect(tc, dir, &t, &a, &b, &c)) {
	  //printf("%f,%f  %f,%f  %f,%f\n", bary.x, a, bary.y, b, bary.z, c);

	  bary = Vec3f(a, b, c);
	  baryImage(x, y, 0) = bary.x;
	  baryImage(x, y, 1) = bary.y;
	  baryImage(x, y, 2) = bary.z;
	}
	
	Vec3f co = mesh.vert[mesh.tris[trii].vi[0]]*bary.x +
	  mesh.vert[mesh.tris[trii].vi[1]]*bary.y +
	  mesh.vert[mesh.tris[trii].vi[2]]*bary.z;
	
	Vec3f no = mesh.norm[mesh.tris[trii].ni[0]]*bary.x +
	  mesh.norm[mesh.tris[trii].ni[1]]*bary.y +
	  mesh.norm[mesh.tris[trii].ni[2]]*bary.z;
	
	uv_co(x, y, 0) = co.x;
	uv_co(x, y, 1) = co.y;
	uv_co(x, y, 2) = co.z;
	
	uv_no(x, y, 0) = no.x;
	uv_no(x, y, 1) = no.y;
	uv_no(x, y, 2) = no.z;
	
	mask(x, y) = 1;

	
	// Create the hi-resolution mesh.
	vertMap(x, y) = submesh.vert.size();
	submesh.restVert.push_back(Vec3f(co.x, co.y, co.z));
	submesh.vert.push_back(Vec3f(co.x, co.y, co.z));
	submesh.tvert.push_back(Vec3f(double(x)/w, double(y)/h, 0.0));
	submesh.vertexSource.push_back(Vec2<int>(x, y));

	// Only add quad-like triangles.
	if (y > 0 && x > 0 && 
	    timage(x-1, y-1) < mesh.tris.size() &&
	    timage(x, y-1) < mesh.tris.size() &&
	    timage(x-1, y) < mesh.tris.size()) {
	  int vi1[3] = {(int)vertMap(x-1, y-1),
			(int)vertMap(x, y-1),
			(int)vertMap(x, y)};
	  int vi2[3] = {(int)vertMap(x-1, y-1),
			(int)vertMap(x, y),
			(int)vertMap(x-1, y)};
	  int ni[3] = {0, 0, 0};
	  submesh.tris.push_back(Mesh::Tri(vi1, ni, vi1));
	  submesh.tris.push_back(Mesh::Tri(vi2, ni, vi2));
	}
      }
    }
  }
   
  return BaseGeometry::ptr(new BaseGeometryType(mesh, uv_co, uv_no, mask, baryImage, timage));

  uv_co = 0;
  uv_no = 0;
  mask = 0;
  
  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      Vec3f o(double(x)/w, double(y)/h, 0);
      Vec3f k(0, 0, 1);

      std::cout << o.x << "," << o.y << "\n";

      for (int k=0; k<(int)triangles.size(); k++) { 
	Vec3d bary;
	double tout;

	if (triangles[k].distance(o) <= 1e-6) {
	  std::cout << "Hit\n";
	}

	if (triangles[k].intersect(o, k, &tout, &bary.x, &bary.y, &bary.z)) {
	  std::cout << "Hit" << x << " " << y << "\n";
	  baryImage(x, y, 0) = bary.x;
	  baryImage(x, y, 1) = bary.y;
	  baryImage(x, y, 2) = bary.z;

	  Vec3f co = mesh.vert[mesh.tris[k].vi[0]]*bary.x +
	    mesh.vert[mesh.tris[k].vi[1]]*bary.y +
	    mesh.vert[mesh.tris[k].vi[2]]*bary.z;
	  
	  Vec3f no = mesh.norm[mesh.tris[k].ni[0]]*bary.x +
	    mesh.norm[mesh.tris[k].ni[1]]*bary.y +
	    mesh.norm[mesh.tris[k].ni[2]]*bary.z;
	  
	  uv_co(x, y, 0) = co.x;
	  uv_co(x, y, 1) = co.y;
	  uv_co(x, y, 2) = co.z;

	  uv_no(x, y, 0) = no.x;
	  uv_no(x, y, 1) = no.y;
	  uv_no(x, y, 2) = no.z;

	  mask(x, y) = 1;
	  break;
	}
      }
    }
  }
  return BaseGeometry::ptr(new BaseGeometryType(mesh, uv_co, uv_no, mask, baryImage, timage));
}



GLuint genDepthTexture(int w,int h){
  GLuint depthTex;
  glGenTextures(1, &depthTex);
  glBindTexture(GL_TEXTURE_2D,depthTex);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,0,0,w,h,0);
  glBindTexture(GL_TEXTURE_2D,0);

  printf("Generated depth texture: %d, %dx%d\n", depthTex,w,h);

  return depthTex;
}


nacb::Imagef getWeightImage(const DisparitySequence & seq, 
			    BaseGeometry::ptr & baseGeometry,
			    DisplacedMesh* mesh,
			    bool occlusion){
  Vec3f cc = seq.E.inverse().submatrix(0, 3, 3, 1);
  Vec3d ccd(cc.x, cc.y, cc.z);

  nacb::Imagef weights(baseGeometry->getWidth(),
		       baseGeometry->getHeight(), 1);
  weights = 1.0;

  // Use the current mesh to get the occlusions.
  if (occlusion && mesh) {
    glViewport(0, 0, seq.image.w, seq.image.h);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    GLuint tex = 0;
    GLuint depthTex = genDepthTexture(seq.image.w, seq.image.h);
    double zNear = 1;
    double zFar = 7.0;

    glGenTextures(1, &tex);

    depthRange(mesh, &(seq.E(2, 0)), zNear, zFar);
    printf("range: %f %f\n", zNear, zFar);
    cameraPerspectiveMatrix(seq.A(0, 0), seq.A(1, 1), 
			    seq.A(0, 2), seq.A(1, 2),
			    seq.image.w, seq.image.h, 
			    zNear, zFar, true);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    nacb::Matrix E = seq.E.transpose();
    glMultMatrixd(E.data);

    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2, 2);

    glColorMask(0,0,0,0);
    glColor3f(1, 1, 1);
    mesh->draw();
    glColorMask(1,1,1,1);

    glDisable(GL_POLYGON_OFFSET_FILL);

    //nacb::Imagef depth(seq.image.w, seq.image.h, 1);
    //glReadPixels(0, 0, depth.w, depth.h, GL_DEPTH_COMPONENT, GL_FLOAT, depth.data);
    //depth.save("/tmp/depth.png");

    glActiveTexture(GL_TEXTURE0);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tex);

    nacb::Image8 image = seq.image.boxFilter(3);
    image.initTexture();
    //glBindTexture(GL_TEXTURE_2D, 0);

    glActiveTexture(GL_TEXTURE1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,depthTex);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    //Grab the depth and put in texture
    glEnable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, depthTex);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,0,0,seq.image.w, seq.image.h, 0);

    glColor3f(1, 1, 1);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_COMPARE_FUNC_ARB, GL_LEQUAL);
    glTexParameteri(GL_TEXTURE_2D,GL_DEPTH_TEXTURE_MODE_ARB,GL_INTENSITY);
    
    glClear(GL_DEPTH_BUFFER_BIT);

    glViewport(0, 0, weights.w, weights.h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();    
    glOrtho(0,1,0,1,-1,1);


    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    for(int i=1;i>=0;i--){
      glActiveTexture(GL_TEXTURE0+i);
      glMatrixMode(GL_TEXTURE);
      glLoadIdentity();
      
      double mat[16]={0.5,0.0,0.0,0.0,
		      0.0,0.5,0.0,0.0,
		      0.0,0.0,0.5,0.0,
		      0.5,0.5,0.5,1.0};
      glMultMatrixd(mat);
      cameraPerspectiveMatrix(seq.A(0, 0), seq.A(1, 1), 
			      seq.A(0, 2), seq.A(1, 2),
			      seq.image.w, seq.image.h, 
			      zNear, zFar, true);
      glMultMatrixd(E.data);
    }

    drawTC(mesh, ccd);

    for(int i=1; i>=0; i--) {
      glActiveTexture(GL_TEXTURE0+i);
      glMatrixMode(GL_TEXTURE);
      glLoadIdentity();
      glDisable(GL_TEXTURE_2D);
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glDisable(GL_TEXTURE_2D);

    glReadPixels(0, 0, weights.w, weights.h, GL_ALPHA, GL_FLOAT, weights.data);
    glDeleteTextures(1, &depthTex);   
    glDeleteTextures(1, &tex);

    // Do a number of smoothing iterations (proportional to size)
    for (int j=seq.image.w; j != 0; j = j >> 1) {
      weights = weights.boxFilter(4); 
    }
    printGLError();
    return weights;
  }

  for (int y=0; y<baseGeometry->getHeight(); y++) {
    for (int x=0; x<baseGeometry->getWidth(); x++) {
      Vec3f no;
      Vec3f co = baseGeometry->displace(x, y, 0, &no);
      Vec3f d = cc - co;
      d.normalize();
      double val = d.dot(no);
      val = std::max(val, 0.0);
      weights(x, y) = (val + weights(x, y))/2;
    }
  }
  weights = weights * baseGeometry->getMask();
  return weights;
}


// Instantiate the functions
template
BaseGeometry::ptr createBaseGeometryFromMesh<ImageBaseGeometry>(Mesh & mesh, int w, int h, 
								DisplacedMesh & submesh);

template
BaseGeometry::ptr createBaseGeometryFromMesh<ImageBaseSkinnedGeometry>(Mesh & mesh, int w, int h, 
								       DisplacedMesh & submesh);
