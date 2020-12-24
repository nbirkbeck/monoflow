#ifndef DEPTH_MESH_H
#define DEPTH_MESH_H

#include <GL/gl.h>
#include <nimage/image.h>
#include <nmath/matrix.h>
#include <boost/function.hpp>


/**
   Abstract representation of a quad-mesh.  Can be 
   used to switch between rendering modes (e.g., display
   lists, vbo's); can be shared also.  
*/
class QuadMesh2D {
 public:
  virtual ~QuadMesh2D(){ }
  virtual void draw() = 0;
  virtual void build() = 0;
  virtual void setSpacing(int sp){
    spacing = sp;
    build();
  }

  virtual int getSpacing() const {
    return spacing;
  }

 protected:
  int spacing;

 private:
  void operator=(QuadMesh2D & other){ }
};


class BasicQuadMesh2D : public QuadMesh2D {
 public:
  BasicQuadMesh2D(int _w, int _h, int _spacing = 2, bool _usingList = false){
    w = _w;
    h = _h;
    spacing = _spacing;
    usingList = _usingList;

    list = 0;
    build();
  }

  ~BasicQuadMesh2D(){
    glDeleteLists(list, 1);
    list = 0;
  }

  virtual void build() {
    if(usingList){
      if(!list)
	list = glGenLists(1);
      glNewList(list, GL_COMPILE);

      draw();
      
      glEndList();
    }
  }

  virtual void draw() {
    int sx = 2, sy = 2;
    for(int y=sy; y<h-1-sy; y+=spacing){
      glBegin(GL_QUAD_STRIP);
      glVertex2f(sx, y);
      glVertex2f(sx, y+spacing);
      
      for(int x=sx+1; x<w-sx; x+=spacing){
	glVertex2f(x, y);
	glVertex2f(x, y+spacing);
      }
      glEnd();
    }
  }
  
 protected:
  int w, h;
  bool usingList;
  GLuint list;
};


class VboQuadMesh2D : public QuadMesh2D {
 protected:
  struct Elements {
    GLenum mode;
    GLint  first;
    GLsizei count;
  };

 public:
  VboQuadMesh2D(int _w, int _h, int _spacing = 2){
    glGenBuffers(1, &vbo);
    w = _w;
    h = _h;
    spacing = _spacing;
    
    build();
  }

  ~VboQuadMesh2D(){
    glDeleteBuffers(1, &vbo);
    vbo = 0;
  }

  void build(){
    int maxv = (int) (ceil(float(h) / spacing) * ceil(float(w) / spacing) * 4);
    nacb::Vec2f * verts = new nacb::Vec2f[maxv];

    elements.clear();

    int sx = 2, sy = 2;
    int nv = 0;

    for(int y=sy; y<h-1-sy; y+=spacing){
      Elements element;
      element.first = nv;
      element.mode = GL_QUAD_STRIP;

      verts[nv++] = nacb::Vec2f(sx, y);
      verts[nv++] = nacb::Vec2f(sx, y + spacing);
      
      for(int x=sx+1; x<w-sx; x+=spacing){
	verts[nv++] = nacb::Vec2f(x, y);
	verts[nv++] = nacb::Vec2f(x, y + spacing);
      }
      element.count = nv - element.first;
      elements.push_back(element);
    }

    glBindBuffer(GL_ARRAY_BUFFER_ARB, vbo);
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, 2*sizeof(float)*nv, verts, GL_DYNAMIC_DRAW_ARB);
    glBindBuffer(GL_ARRAY_BUFFER_ARB, 0);

    delete [] verts;
  }

  virtual void draw() {
    glBindBuffer(GL_ARRAY_BUFFER_ARB, vbo);
    glEnableClientState(GL_VERTEX_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, 0);

    for(int i=0; i<elements.size(); i++)
      glDrawArrays(elements[i].mode, elements[i].first, elements[i].count);
    
    glDisableClientState(GL_VERTEX_ARRAY);
  }
  
 protected:

  std::vector<Elements> elements;
  GLuint vbo;
  int w, h;
};



class DepthMesh{
public:
  DepthMesh(const nacb::Image8 & image,
	    const nacb::Imagef & depth,
	    const nacb::Matrix & P);

  ~DepthMesh();

  void draw(GLuint texUse = 0, GLuint depthTexUse = 0);
  void setShading(float s);
  float getShading();
  void setSpacing(int s);
  int getSpacing();

  void setDepth(const nacb::Imagef & depth);

protected:
  void setProgramParameters();
  void multModelviewMatrix();
  
  GLuint depthTex, tex;
  
  nacb::Matrix E, K;
  int spacing;
  
  int w; int h;

  ///The programs could be static.
  GLuint prog, fragShader, vertShader;

  float shading;

  QuadMesh2D * mesh;
};


// The following functions are used to create mesh representations
// similar to a depth mesh (but not necessarily a depth mesh).
// They are in this file, as it is currently the most relevant.

double fromWindowDepth(double zwin, double minz, double maxz);

nacb::Imagef renderDepth(const nacb::Matrix & K, const nacb::Matrix & E,
			 int w, int h, double nr, double fr,
			 const boost::function<void()> & drawfunc);

void renderDepthMesh(const nacb::Matrix & K, const nacb::Matrix & E,
		     int w, int h,
		     double nr, double fr,
		     const boost::function<void()> & drawFunc, 
		     std::vector<nacb::Vec3<int> > & tris, 
		     std::vector<nacb::Vec3f> & vert);

#endif //DEPTH_MESH_H
