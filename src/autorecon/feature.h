#ifndef FEATURE_H
#define FEATURE_H

#include <nimage/image.h>
#include <nmath/matrix.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include <vector>
#include <stdio.h>

using nacb::Vec3f;
using nacb::Vec2f;
using namespace std;
using nacb::Imagef;
using nacb::Image8;
using nacb::Matrix;

enum{BASIC_FEATURE,IMAGE_FEATURE,SIFT_FEATURE};

class Feature{
public:
  unsigned char type;
  Vec2f point;
  Feature(unsigned char type=0,const Vec2f & point=Vec2f(0,0)){
    this->type=type;
    this->point=point;
  }
  virtual ~Feature(){
  }
  virtual double compare(Feature * feature){
    return 10e10;
  }
  static Feature * loadFeature(void * buffer,int & len);
  static Feature * loadFeature(FILE * file);
  virtual int load(void * buffer,int len){
    if(len>=8){
      memcpy(point.data,buffer,sizeof(float)*2);
      return 8;
    }
    return 0;
  }
  virtual int save(void * buffer){
    if(buffer==0)return sizeof(float)*2+1;
    unsigned char * cbuffer=(unsigned char *)buffer;
    cbuffer[0]=type;
    memcpy(cbuffer+1,point.data,sizeof(float)*2);
    return sizeof(float)*2+1;
  }
  virtual int load(FILE * file){
    return 2==fread(point.data,sizeof(float),2,file);
  }
  virtual int save(FILE * file){
    if(!fwrite(&type,sizeof(unsigned char),1,file))return 0;
    return 2==fwrite(point.data,sizeof(float),2,file);
  }
};

class ImageFeature: public Feature{
public:
  static int correlation;
  nacb::Imagef image;
  ImageFeature(const Vec2f & point=Vec2f(0,0),
	       const nacb::Imagef & image=nacb::Imagef(1,1,1)):Feature(IMAGE_FEATURE,point){
    this->image=image;
    
    double mn=0;
    for(int y=0;y<image.h;y++)
      for(int x=0;x<image.w;x++)
	for(int k=0;k<image.nchannels;k++)
	  mn+=image(x,y,k);
    
    mn/=(image.w*image.h*image.nchannels);
    double len=0;
    for(int y=0;y<image.h;y++){
      for(int x=0;x<image.w;x++){
	for(int k=0;k<image.nchannels;k++){
	  image(x,y,k)-=mn;
	  double v=image(x,y,k);
	  len+=v*v;
	}
      }
    }
    if(ImageFeature::correlation){
      len=sqrt(len);
      //remove the length
      for(int y=0;y<image.h;y++){
	for(int x=0;x<image.w;x++){
	  for(int k=0;k<image.nchannels;k++){
	    image(x,y,k)/=len;
	  }
	}
      } 
    }
  }
  virtual ~ImageFeature(){}

  virtual double compare(Feature * feature){
    if(feature->type==IMAGE_FEATURE){
      ImageFeature * other=(ImageFeature *)feature;
      double sc=0;
      if(ImageFeature::correlation){
	double sc=0;
	for(int y=0;y<image.h;y++){
	  for(int x=0;x<image.w;x++){
	    for(int k=0;k<image.nchannels;k++){
	      sc+=image(x,y,k)*other->image(x,y,k);
	    }
	  }
	}
	return -sc;
      }
      else{
	for(int y=0;y<image.h;y++){
	  for(int x=0;x<image.w;x++){
	    for(int k=0;k<image.nchannels;k++){
	      double d=image(x,y,k)-other->image(x,y,k);
	      sc+=d*d;
	    }
	  }
	}
	return sc;
      }
    }
    else return Feature::compare(feature);
  }
  virtual int load(void * buffer,int len){
    int parent_len;
    if(!(parent_len=Feature::load(buffer,len)))return 0;
    len-=parent_len;

    if(len<12)return 0;
    
    int dims[3];
    unsigned char * cbuffer=(unsigned char *)buffer;
    memcpy(dims,cbuffer+parent_len,sizeof(int)*3);
    int required=dims[0]*dims[1]*dims[2]*4;


    len-=sizeof(int)*3;
    if(len<required)return 0;
    
    image=nacb::Imagef(dims[0],dims[1],dims[2]);
    memcpy(image.data,cbuffer+parent_len+sizeof(int)*3,required);

    return parent_len+sizeof(int)*3+required;
  }
  virtual int save(void * buffer){
    int parent_len=Feature::save(buffer);
    int len=parent_len+3*sizeof(int)+image.w*image.h*image.nchannels*sizeof(float);
    if(buffer==0)return len;
    
    //need single bytes
    unsigned char * cbuffer=(unsigned char *)buffer;
    int dims[3]={image.w,image.h,image.nchannels};
    memcpy(cbuffer+parent_len,dims,sizeof(int)*3);
    memcpy(cbuffer+parent_len+sizeof(int)*3,image.data,
	   sizeof(float)*image.w*image.h*image.nchannels);
    return len;
  }
  virtual int load(FILE * file){
    if(!Feature::load(file))return 0;
    int dims[3];
    if(3!=fread(dims,sizeof(int),3,file))return 0;
    image=nacb::Imagef(dims[0],dims[1],dims[2]);
    if(dims[0]*dims[1]*dims[2]!=
       (int)fread(image.data,sizeof(float),dims[0]*dims[1]*dims[2],file))
      return 0;
    return 1;
  }
  virtual int save(FILE * file){
    if(!Feature::save(file))return 0;
    int dims[3]={image.w,image.h,image.nchannels};
    if(3!=fwrite(dims,sizeof(int),3,file))return 0;
    if(dims[0]*dims[1]*dims[2]!=
       (int)fwrite(image.data,sizeof(float),dims[0]*dims[1]*dims[2],file))
      return 0;
    return 1;
  }
};


std::vector<Feature*> getFeatures(nacb::Image8 & image);

//Get matching feature points (matrix form) for all matches, returns empty vector if no-matches.
std::vector<nacb::Matrix> getMatchedPoints(const std::vector<Feature *> & f1,
					   const std::vector<Feature *> & f2,
					   const std::vector<int> & matches);


#endif
