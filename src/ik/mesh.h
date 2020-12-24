#ifndef MESH_H
#define MESH_H

#ifdef USING_GLEW
#include <GL/glew.h>
#endif

#include <vector>
#include <set>
#include <map>
#include <GL/gl.h>
#include <nmath/matrix.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include <assert.h>
#include <nmath/mat4x4.h>
#include <nimage/image.h>
using namespace std;

class Mesh{
public:
  class Face{
  public:
    int vi[4];
    int ni[4];
    int tci[4];
    Face(){
      memset(this,0,sizeof(Face));
    }
    Face(int vi[],int ni[]){
      memset(this,0,sizeof(Face));
      memcpy(this->vi,vi,sizeof(this->vi));
      memcpy(this->ni,ni,sizeof(this->ni));
    }
    Face(int vi[],int ni[],int tci[]){
      memcpy(this->vi,vi,sizeof(this->vi));
      memcpy(this->ni,ni,sizeof(this->ni));
      memcpy(this->tci,tci,sizeof(this->tci));
    }
    Face(int v0,int n0,
	 int v1,int n1,
	 int v2,int n2,
	 int v3,int n3){
      memset(this,0,sizeof(Face));

      vi[0]=v0;
      vi[1]=v1;
      vi[2]=v2;
      vi[3]=v3;
      
      ni[0]=n0;
      ni[1]=n1;
      ni[2]=n2;
      ni[3]=n3;
    }
    Face(int v0,int n0,int tc0,
	 int v1,int n1,int tc1,
	 int v2,int n2,int tc2,
	 int v3,int n3,int tc3){
      vi[0]=v0;
      vi[1]=v1;
      vi[2]=v2;
      vi[3]=v3;
      
      ni[0]=n0;
      ni[1]=n1;
      ni[2]=n2;
      ni[3]=n3;

      tci[0]=tc0;
      tci[1]=tc1;
      tci[2]=tc2;
      tci[3]=tc3;
    }
    bool operator==(const Face & other)const{
      return memcmp(this,&other,sizeof(Face));
    }
  };
  class Tri{
  public:
    int vi[3];
    int ni[3];
    int tci[3];
    Tri(){
      memset(this,0,sizeof(Tri));
    }
    Tri(int vi[],int ni[]){
      memset(this,0,sizeof(Tri));
      memcpy(this->vi,vi,sizeof(this->vi));
      memcpy(this->ni,ni,sizeof(this->ni));
    }
    Tri(int vi[],int ni[],int tci[]){
      memcpy(this->vi,vi,sizeof(this->vi));
      memcpy(this->ni,ni,sizeof(this->ni));
      memcpy(this->tci,tci,sizeof(this->tci));
    }
    bool operator==(const Tri & other)const{
      return memcmp(this,&other,sizeof(Tri));
    }
  };

  struct bone_weight_t{
    float weight;
    short bone;
    bone_weight_t(short b=0,float w=0):weight(w),bone(b){}

    bool operator==(const bone_weight_t & o)const{
      return weight==o.weight && bone==o.bone;
    }
  };
  typedef std::vector<bone_weight_t>  weight_vector_t;
  typedef std::vector<weight_vector_t> bone_weights_t;

  struct bounding_sphere_t{
    nacb::Vec3f c;
    nacb::Vec3f restc;
    double r;
    vector<bone_weight_t> weights;
    vector<short> interacts;

    bounding_sphere_t(const nacb::Vec3f & center=nacb::Vec3f(0,0,0),
		     double radius=0,
		      const vector<bone_weight_t> & w=vector<bone_weight_t>(),
		      const vector<short> & inters=vector<short>())
    : c(center),  restc(center), r(radius), weights(w), interacts(inters){ 
    }
    void cleanInteracts(){
      std::sort(interacts.begin(), interacts.end());
      vector<short>::iterator newend=unique(interacts.begin(), interacts.end());
      printf("copying results\n");
      vector<short> newinter(newend-interacts.begin());

      std::copy(interacts.begin(), newend, newinter.begin());
      interacts=newinter;
    }
    bool operator==(const bounding_sphere_t & other)const{
      return false;
    }
  };


  vector<nacb::Vec3f> tvert;
  vector<nacb::Vec3f> vert;
  vector<nacb::Vec3f> norm;
  vector<nacb::Vec3f> restVert;
  
  vector<Face> poly;
  vector<Tri> tris;

  bone_weights_t  bone_weights;

  vector<bounding_sphere_t> spheres;
  float * vinds;
  int vinds_size;
  unsigned int vbo;

  Mesh(){
    vinds = 0;
    vinds_size = 0;
    vbo = 0;
  }

  Mesh(const Mesh & m){
    vinds= 0;
    vinds_size = 0;
    vbo = 0;
    operator=(m);
  }

  Mesh(const nacb::Imagef & coords, const nacb::Image8 & mask);

  ~Mesh(){
    delete [] vinds;
    vinds_size = 0;
    vbo = 0;
  }
  void initNormalsFlat();
  void initNormals();
  void initNormalsAutoSmooth(double degreeThresh);

  int loadobj(const char * fname,vector<string> * unused=0);
  int saveobj(const char * fname,bool append=false, bool deformed=false);
  
  void scale(const nacb::Vec3f & sc){
    for(int i=0; i<(int)restVert.size(); i++){
      restVert[i].x*=sc.x;
      restVert[i].y*=sc.y;
      restVert[i].z*=sc.z;
    }
    vert=restVert;
  }
  void translate(const nacb::Vec3f & tr){
    for(int i=0; i<(int)restVert.size(); i++){
      restVert[i].x+=tr.x;
      restVert[i].y+=tr.y;
      restVert[i].z+=tr.z;
    }
    vert=restVert;
  }

  template <class T> void animate(T);

  void draw(bool fancy=false);

  void drawWithTC(const std::vector<nacb::Vec3f> & tc);

  //Should be called before drawVertexElements and after any vertex changes.
  void initVertexElements();
  void drawVertexElements();
  bool hasGeometry(){
    return vert.size()>0;
  }

  void resetSpheres(){
    spheres.clear();
  }
  void addSphere(const nacb::Vec3f & p,
		 const double r,
		 const vector<bone_weight_t> & weights,
		 const vector<short> & interacts=vector<short>()){
    spheres.push_back(bounding_sphere_t(p,r,weights,interacts));
  }
  Mesh subdivide() const;
  
  Mesh keepVert(const vector<int> & vkeep) const;

  const Mesh & operator=(const Mesh & mesh){
    if(this == &mesh)return *this;

    tvert = mesh.tvert;
    vert = mesh.vert;
    norm = mesh.norm;
    restVert = mesh.restVert;
    poly = mesh.poly;
    tris = mesh.tris;
    bone_weights = mesh.bone_weights;
    spheres = mesh.spheres;
    if(vinds) delete [] vinds;
    vinds = 0;
    vinds_size = 0;
    return *this;
  }
  //Only adds faces
  void fillSimpleHoles();

  std::set<int> findBorderVertices() const {
    // Find any edges without opposites.
    std::set<int> borderVerts;
    
    std::map<int, int> edgeCount;
    
    for (int i=0; i<(int)tris.size(); i++) {
      for (int k=0; k<3; k++) {
	int v1 = tris[i].vi[k];
	int v2 = tris[i].vi[(k + 1)%3];	
	int index = std::max(v1, v2)*vert.size() + std::min(v1, v2);
	
	if (!edgeCount.count(index)) {
	  edgeCount[index]++;
	}
	else
	  edgeCount[index] = 1;
      }
    }

    std::map<int, int>::iterator it = edgeCount.begin();
    std::map<int, int>::iterator last = edgeCount.end();
    
    while (it != last) {
      if (it->second == 1) {
	borderVerts.insert(it->first / vert.size());
	borderVerts.insert(it->first % vert.size());
      }
      it++;
    }
    return borderVerts;
  }
};

bool operator==(const vector<Mesh::bone_weight_t> & v1,
		const vector<Mesh::bone_weight_t> & v2);

void intersectRayWithMesh(const Mesh & mesh,
			  const nacb::Vec3f & o,
			  const nacb::Vec3f & k,
			  double & mindist,
			  double & maxdist);

void intersectRayWithMesh(const Mesh & mesh,
			  const nacb::Vec3f & o, const nacb::Vec3f & k,
			  double & mindist,
			  double & maxdist,
			  int & triangleIndex,
			  nacb::Vec3f & bary, 
			  bool returnNear = true);


//Used in warp texture
void drawTC(Mesh * mesh, const nacb::Vec3d & ccd);

//Used in warp texture
void depthRange(Mesh * mesh, const double * zdir, 
		double & zmin, double & zmax);



/**
   Return the bary-centric coordinate image for the mesh 
   of size tw, th (zeros where there is no triangle).

   If useTexCoords, then it is for texture space, otherwise
   the projection matrices must be set and it is for vertices.

   If timage!=null then fill it with the triangle indices (-1 
   where there is no triangle).
*/
nacb::Imagef getBaryImage(Mesh & mesh, int tw, int th, 
			  bool useTexCoords = true, nacb::Image32 * timage = 0);

#endif
