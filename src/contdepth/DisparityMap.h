#ifndef DISPARITY_MAP_H
#define DISPARITY_MAP_H

#include <nimage/image.h>
#include <nmath/matrix.h>
#include <assert.h>


/**
   \brief A class to maintain a discrete disparity map (that
          maps to corresponding inverse depth.
*/
template <class T>
class DisparityMap {

 public:
  DisparityMap(int w, 
	       int h,
	       double _dmin, 
	       double _dmax, int _nd) : dmin(_dmin), dmax(_dmax), nd(_nd){
    labels = nacb::Image<T>(w, h, 1);        
  }

  template <class T2>
    DisparityMap(DisparityMap<T2> & dmap) : dmin(dmap.dmin), dmax(dmap.dmax), nd(dmap.nd) {
    
    labels = nacb::Image<T>(dmap.labels.w, dmap.labels.h, 1);

    for(int y=0; y<labels.h; y++){
      for(int x=0; x<labels.w; x++){
	labels(x, y) = (T)dmap.labels(x, y);
      }
    }
  }

  template <class T2>
  float labelToDisparity(const T2 & label) const {
    return float(label)*(dmax - dmin)/float(nd - 1) + dmin;
  }

  float labelToDisparitySlope() const {
    return (dmax - dmin)/float(nd - 1);
  }
  
  T disparityToLabel(const float & d) const {
    int index = round((nd-1)*(d - dmin)/dmax);
    return std::max(0, std::min(index, nd - 1));
  }
  
  const float bilinear(float x, float y) const {
    if(x < 0 || y < 0 || x >= labels.w || y >= labels.h)
      return -1;

    int x0 = int(x), y0 = int(y);
    int x1 = std::min(labels.w - 1, x0 + 1);
    int y1 = std::min(labels.h - 1, y0 + 1);
    float a = x - x0, b = y - y0;

    return float(getDisparity(x0, y0))*(1.0f - a)*(1.0f - b) +
      float(getDisparity(x1, y0))*a*(1.0f - b) +
      float(getDisparity(x1, y1))*a*b +
      float(getDisparity(x0, y1))*(1.0f - a)*b;      
  }

  DisparityMap& operator=(const nacb::Matrix & matrix){
    assert(matrix.m == labels.width*labels.height);
    
    for(int i=0; i<matrix.m; i++)
      labels.data[i] = matrix.data[i];

    return *this;
  }

  const T & operator[](int i) const {
    return labels.data[i];
  }

  T & operator[](int i) {
    return labels.data[i];
  }

  const T & operator()(int x, int y) const {
    return labels(x, y);
  }

  T & operator()(int x, int y) {
    return labels(x, y);
  }

  float getDisparity(int x, int y) const {
    return labelToDisparity(labels(x, y));
  }

  nacb::Imagef getDepth() const {
    assert(0);
    return nacb::Imagef();
  }

  nacb::Imagef getDepth(const nacb::Image8 & mask, const T & defaultLabel) const {
    assert(mask.w == labels.w && mask.h == labels.h);
    nacb::Imagef depth(labels.w, labels.h, 1);

    for(int y=0; y<labels.h; y++){
      for(int x=0; x<labels.w; x++){
	T label = defaultLabel;
	if(mask(x, y))
	  label = labels(x, y);

	depth(x, y) = 1.0/labelToDisparity(label);
      }
    }
    return depth;
  }

  nacb::Imagef getMappedDisparity() const {
    nacb::Imagef disp(labels.w, labels.h, 1);
    for(int y=0; y<labels.h; y++){
      for(int x=0; x<labels.w; x++){
	disp(x, y) = labelToDisparity(labels(x, y));
      }
    }
    return disp;
  }

  nacb::Imagef getNormalizedDisparityImage() const {
    nacb::Imagef l(labels.width, labels.height, 1);
    for(int y=0; y<l.height; y++){
      for(int x=0; x<l.width; x++){
	l(x, y) = labels(x, y);
      }
    }
    l *= (1.0/(nd - 1));
    return l;
  }

  nacb::Matrix getFlatMatrix() const {
    nacb::Matrix m(labels.width*labels.height, 1);
    for(int i=0; i<labels.width*labels.height; i++)
      m.data[i] = labels.data[i];
    return m;
  }

  void setFromDisparity(nacb::Imagef & disp) {
    if(disp.w != labels.w ||
       disp.h != labels.h)
      throw "incompatible image dimensions.";
    
    for(int y=0; y<disp.h; y++){
      for(int x=0; x<disp.w; x++){
	labels(x, y) = (T)disparityToLabel(disp(x, y));
      }
    }
  }

  nacb::Image<T> labels;

  float dmin, dmax;  
  int  nd;
};

typedef DisparityMap<uint16_t> DiscreteDisparityMapInt;

#endif // DISPARITY_MAP_H
