#ifndef EXPLICIT_DEPTH_OFFSET_MESH_H
#define EXPLICIT_DEPTH_OFFSET_MESH_H

#include <GL/glu.h>

#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include <nimage/image.h>
#include "utigl/glcommon.h"

#include <vector>

#include "autorecon/recon_globals.h"
#include "autorecon/recon_geometry.h"

using namespace nacb;
using namespace std;


/**
   A class for rendering from depth and offset meshes 
   at different stages of its displacement.   
*/
class ExplicitDepthOffsetMesh {
public:
  nacb::Image8 image;
  nacb::Imagef depth;
  nacb::Imagef  flow;
  
  Vec3f colorMin;
  Vec3f colorMax;

  vector<Vec3f> coords;
  double scaleOffs;
  Matrix E, K;
  bool visible;
  GLuint tex;
  
  Vec3f cc;

  ExplicitDepthOffsetMesh() { }

  ExplicitDepthOffsetMesh(const Matrix & K_,
			  const Matrix & E_,
			  const nacb::Image8 & image_,
			  const nacb::Imagef & depth_,
			  const nacb::Imagef & flow_) {
    E = E_;
    K = K_;

    image = image_;
    depth = depth_;
    flow = flow_;

    visible = true;
    
    scaleOffs = 1;
    colorMin = Vec3f(-1, -1, -1)*0.01;
    colorMax = Vec3f( 1,  1,  1)*0.01;
    
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    image.initTexture();
    glBindTexture(GL_TEXTURE_2D, 0);

    cc = E.inverse().submatrix(0, 3, 3, 1);

    computeCoords();
  }

  void computeCoords(){
    Matrix P = K*Matrix::eye(3, 4)*E;
    Matrix P3x3inv = P.submatrix(0, 0, 3, 3).inverse();
    Matrix KRt = P.submatrix(0, 3, 3, 1);
    
    coords = std::vector<Vec3f>(image.w*image.h);

    for(int y=0; y<image.h-1; y++){
      for(int x=0; x<image.w-1; x++){
	Vec3f co = backProject(P3x3inv, KRt, (float)x, (float)y, depth(x, y));
	coords[y*image.width + x] = co;
      }
    }
  }

  Vec3f getCoordsOffset(int x, int y) {
    Vec3f offs(flow(x, y, 0), flow(x, y, 1), flow(x, y, 2));
    return coords[y*image.width + x] + offs*scaleOffs;
  }

  Vec3f getCoords(int x, int y) {
    return coords[y*image.width + x];
  }
 
  Vec3f getFlowColor(int x, int y) const {
    Vec3f color(flow(x, y, 0), flow(x, y, 1), flow(x, y, 2));
    color = (color - colorMin);
    Vec3f range = (colorMax - colorMin) ;
    color.x *= (1.0/range.x);
    color.y *= (1.0/range.y);
    color.z *= (1.0/range.z);    
    
    return color;
  }

  void drawTriangles(int skip, Vec3f (ExplicitDepthOffsetMesh::* getCoordsFunc)(int x, int y), 
		     bool colors = false){
    glBegin(GL_TRIANGLES);
    
    for(int y=0; y<image.h-skip; y+=skip){
      for(int x=0; x<image.w-skip; x+=skip){
	if(image(x, y, 3) >= 100 &&
	   image(x+skip, y, 3) >= 100 &&
	   image(x, y+skip, 3) >= 100 &&
	   image(x+skip, y+skip, 3) >= 100){
	  
	  Vec3f e1 = (this->*getCoordsFunc)(x+skip, y) - (this->*getCoordsFunc)(x, y);
	  Vec3f e2 = (this->*getCoordsFunc)(x+skip, y+skip) - (this->*getCoordsFunc)(x, y);
	  Vec3f n = e1 ^ e2;

	  n.normalize();

	  Vec3f diff = cc - (this->*getCoordsFunc)(x+skip, y);
	  diff.normalize();

	  if(diff.dot(n) <= -0.1){
	    if(colors)
	      glColor3fv(getFlowColor(x, y).data);
	    
	    glNormal3fv(n.data);
	    glTexCoord3fv(getCoords(x, y).data);
	    glVertex3fv((this->*getCoordsFunc)(x, y).data);
	    
	    if(colors)
	      glColor3fv(getFlowColor(x + skip, y).data);
	    
	    glTexCoord3fv(getCoords(x+skip, y).data);
	    glVertex3fv((this->*getCoordsFunc)(x+skip, y).data);
	    
	    if(colors)
	      glColor3fv(getFlowColor(x + skip, y + skip).data);
	    
	    glTexCoord3fv(getCoords(x+skip, y+skip).data);
	    glVertex3fv((this->*getCoordsFunc)(x+skip, y+skip).data);
	  }
	  e1 = (this->*getCoordsFunc)(x+skip, y+skip) - (this->*getCoordsFunc)(x, y);
	  e2 = (this->*getCoordsFunc)(x, y+skip) - (this->*getCoordsFunc)(x, y);
	  n = e1 ^ e2;

	  n.normalize();

	  if(diff.dot(n) <= -0.1){
	  
	    if(colors)
	      glColor3fv(getFlowColor(x, y).data);
	  
	    glTexCoord3fv(getCoords(x, y).data);
	    glVertex3fv((this->*getCoordsFunc)(x, y).data);

	    if(colors)
	      glColor3fv(getFlowColor(x+skip, y+skip).data);

	    glTexCoord3fv(getCoords(x+skip, y+skip).data);
	    glVertex3fv((this->*getCoordsFunc)(x+skip, y+skip).data);

	    if(colors)
	      glColor3fv(getFlowColor(x, y+skip).data);

	    glTexCoord3fv(getCoords(x, y+skip).data);
	    glVertex3fv((this->*getCoordsFunc)(x, y+skip).data);
	  }
	}
      }
    }
    glEnd();
  }

  void drawColors(int skip){
    drawTriangles(skip, &ExplicitDepthOffsetMesh::getCoords, true);
  }

  void draw(int skip = 4, bool showOffs = true, bool showColors = false){
    if(!visible)return;

    if(showColors){
      drawColors(skip);
      return;
    }

    Vec3f (ExplicitDepthOffsetMesh::* getCoordsFunc)(int x, int y) = 
      &ExplicitDepthOffsetMesh::getCoordsOffset;
    
    if(!showOffs)
      getCoordsFunc = &ExplicitDepthOffsetMesh::getCoords;

    glBindTexture(GL_TEXTURE_2D, tex);

    glMatrixMode(GL_TEXTURE);
    glPushMatrix();
    glLoadIdentity();
    glScalef(0.5, 0.5, 0.5);
    glTranslatef(1, 1, 1);
    cameraPerspectiveMatrix(K(0,0), K(1,1), K(0,2), K(1,2), image.w, image.h, 0.1, 100, 1);
    glMultMatrixd(E.transpose().data);

    drawTriangles(skip, getCoordsFunc, false);
    
    glMatrixMode(GL_MODELVIEW);


    if(!showOffs){
      GLUquadric * cylinder = gluNewQuadric();
      GLUquadric * cone = gluNewQuadric();

      glPushAttrib(GL_ALL_ATTRIB_BITS);
      glDisable(GL_TEXTURE_2D);
      glEnable(GL_LIGHTING);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_NORMALIZE);
      glDepthMask(0);
      skip = std::max(8, skip*4);

      for(int its=0; its<2; its++){
	if (its == 1)
	  glDepthFunc(GL_GREATER);
	else
	  glDepthFunc(GL_LESS);

	glColor3f(its==1, its==0, 0);
	
	//glBegin(GL_LINES);
	for(int y=0; y<image.h-skip; y+=skip){
	  for(int x=0; x<image.w-skip; x+=skip){
	    if(image(x, y, 3) >= 100 &&
	       image(x+skip, y, 3) >= 100 &&
	       image(x, y+skip, 3) >= 100 &&
	       image(x+skip, y+skip, 3) >= 100){
	  
	   	    
	      Vec3f offs = Vec3f(flow(x, y, 0), flow(x, y, 1), flow(x, y, 2));
	      offs *= scaleOffs;
	      Vec3f n = offs;

	      double len = n.normalize();
	      Vec3f co = getCoords(x, y);

	      float colen = co.len();
	      if(!isnormal(colen))
		continue;
	      

	      glPushMatrix();
	      glTranslatef(co.x, co.y, co.z);
	      glScalef(len, len, len);

	      glRotatef(atan2(n.y, n.x)*180.0/M_PI, 0, 0, 1);
	      glRotatef(acos(n.z)*180.0/M_PI, 0, 1, 0);
	
	      gluCylinder(cylinder, 0.1, 0.1, 0.75, 10, 10);
	      glTranslatef(0, 0, 0.75);
	      gluCylinder(cone, 0.15, 0.0, 0.25, 10, 10);
	      glPopMatrix();
	      /*
		n.normalize();
		//glNormal3fv(n.data);
	    
		glVertex3fv(getCoords(x, y).data);	    
		offs *= scaleOffs;
		offs += getCoords(x, y);
		glVertex3fv(offs.data);
	      */
	    }
	  }
	}
      }
      //glEnd();

      gluDeleteQuadric(cylinder);
      gluDeleteQuadric(cone);

      glDepthFunc(GL_LESS);

      glPopAttrib(); 
    }

    glMatrixMode(GL_TEXTURE);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
  }

  
  static std::vector<ExplicitDepthOffsetMesh> load(const char * imageName,
						   const char * calibName,
						   const char * depthName,
						   const char * flowName, int n = 8) {
    std::vector<ExplicitDepthOffsetMesh> meshes;

    for(int k=0; k<n; k++) {
      ExplicitDepthOffsetMesh mesh;
      Matrix K = Matrix(3, 3);
      Matrix E = Matrix(4, 4);
      Matrix P = Matrix(3, 4), d = Matrix(5, 1);

      try {
	ClbFile::read(((boost::format(calibName) % k).str().c_str()), K.data, E.data, P.data, d.data);
      
	nacb::Image8 image = nacb::Image8((boost::format(imageName) % k).str().c_str());
	nacb::Imagef depth = nacb::Imagef((boost::format(depthName) % k).str().c_str());
	nacb::Imagef flow(depth.w, depth.h, 3);

	if(flowName)
	  flow = nacb::Imagef((boost::format(flowName) % k).str().c_str());
	else
	  flow = 0;

	mesh = ExplicitDepthOffsetMesh(K, E, image, depth, flow);
      }
      catch (exception & e) {
	std::cout << e.what();
	
	ClbFile::read(calibName, K.data, E.data, P.data, d.data);

	nacb::Image8 image = nacb::Image8(imageName);
	nacb::Imagef depth = nacb::Imagef(depthName);
	nacb::Imagef flow(depth.w, depth.h, 3);

	if(flowName)
	  flow = nacb::Imagef(flowName);
	else
	  flow = 0;

	mesh = ExplicitDepthOffsetMesh(K, E, image, depth, flow);
      }
      
      mesh.scaleOffs = 1;

      meshes.push_back(mesh);
    }
    return meshes;
  }


  static void drawBlendedMeshes(std::vector<ExplicitDepthOffsetMesh> & meshes,
				double t, int skip, bool showOffs, bool backward, bool drawColors) {
    int lo = (int) floor(t);
    int hi = (int) ceil(t);
    
    // Use the lower of the two if doing forward-backward, otherwise round.
    int ind = backward ? int(lo) : (int)round(t);
    double a = t - ind;
    glColor4f(1, 1, 1, 1.0);//1.0 - a);
    
    meshes[ind].scaleOffs = a;
    meshes[ind].draw(skip, showOffs, drawColors);
    
    
    if (backward && hi < (int)meshes.size()) {
      glClear(GL_DEPTH_BUFFER_BIT);
      
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      
      glColor4f(1, 1, 1, a);
      
      meshes[hi].scaleOffs = -(1.0 - a);
      meshes[hi].draw(skip, showOffs, drawColors);
      
      glColor3f(1, 1, 1);
    }
    
  }

};


#endif //EXPLICIT_DEPTH_OFFSET_MESH
