#ifndef FLOW_MAP_H
#define FLOW_MAP_H

#include <nmath/vec3.h>

class FlowMap {

 public:
  FlowMap(){

  }
};


class DiscreteFlowMap : public FlowMap {
  
 public:
  DiscreteFlowMap(int w, int h, float _dmaxx, float _dmaxy, float _dmaxz,
		  int _nx, int _ny, int _nz): dmaxx(_dmaxx), dmaxy(_dmaxy), dmaxz(_dmaxz), 
                                              nx(_nx), ny(_ny), nz(_nz), labels(w, h, 1)
                                              
  {
    swinx = nx*2 + 1;
    swiny = ny*2 + 1;
    swinz = nz*2 + 1;
  }

  nacb::Vec3f getOffset(int x, int y) const {
    return getOffsetForLabel(labels(x, y));
  }

  nacb::Vec3f getOffsetForLabel(int i) const {
    nacb::Vec3f di((i % swinx) - nx, 
		   (i/swinx)%swiny - ny, 
		   i/(swinx*swiny) - nz);
    
    di.x *= (dmaxx/std::max(nx, 1));
    di.y *= (dmaxy/std::max(ny, 1));
    di.z *= (dmaxz/std::max(nz, 1));

    return di;
  }
  
  const unsigned int & operator()(int x, int y) const {
    return labels(x, y);
  }

  unsigned int & operator()(int x, int y){
    return labels(x, y);
  }

  int getLabel(int x, int y) const {
    return labels(x, y);
  }
  
  int getNumLabels() const {
    return swinx*swiny*swinz;
  }

  int getZeroLabel() const {
    return nz*swinx*swiny + ny*swinx + nx;
  }

 protected:
  float dmaxx, dmaxy, dmaxz;
  int nx, ny, nz;
  int swinx, swiny, swinz;
  nacb::Image<unsigned int> labels;
};



#endif // FLOW_MAP_H
