/**
   These used to be scattered about.  Mainly image processing things now.
*/
#ifndef IK_FILTERS_H
#define IK_FILTERS_H

#include <nimage/image.h>
#include <nmath/vec3.h>

nacb::Image8 bilateralSilEdges(nacb::Image8 & image, int hwin);


template <class T>
nacb::Imagef grad(nacb::Image<T> & im,int chan=0);

/*
  \brief Compute motion flow vectors between the two images.
  \param im1 The first image.
  \param im2 The second image.
  \param hwin The half window size for neighborhood
 */
nacb::Imagef motion_flow(nacb::Image8 & im1,nacb::Image8 & im2,int hwin=3,int its=3);

template <class T>
nacb::Image<T> box_filt(nacb::Image<T> & im,int hs);


template<class T>
void bilinear(const nacb::Image<T> & im, float x, float y, nacb::Vec3f & col){
  int x0=std::min(std::max((int)x,0),im.w-1);
  int y0=std::min(std::max((int)y,0),im.h-1);
  int x1=std::min(x0+1,im.w-1);
  int y1=std::min(y0+1,im.h-1);
  double a=x-x0;
  double b=y-y0;
 
  int k=std::min(3,im.nchannels);
  for(int i=0;i<k;i++){
    col.data[i]=(1.0-a)*(1.0-b)*((float)im(x0,y0,i))+
      (a)*(1.0-b)*((float)im(x1,y0,i))+
      (a)*(b)*((float)im(x1,y1,i))+
      (1.0-a)*(b)*((float)im(x0,y1,i));
  }
}

#endif
