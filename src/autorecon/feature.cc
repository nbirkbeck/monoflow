#include <assert.h>
#include "feature.h"
#include "sift.h"
#include "recon_globals.h"

int ImageFeature::correlation=0;

Feature * Feature::loadFeature(void * buffer,int & len){
  unsigned char type;
  unsigned char * cbuffer=(unsigned char*)buffer;
  type=cbuffer[0];
  cbuffer++;

  Feature * feat;
  if(type==BASIC_FEATURE)feat=new Feature(type);
  else if(type==IMAGE_FEATURE)feat=new ImageFeature();
  else if(type==SIFT_FEATURE)feat=new SiftFeature();
  else {
    len=0;
    return 0;
  }
  
  if(len=feat->load(cbuffer,len))return feat;

  len=0;
  delete feat;
  return 0;
}

Feature * Feature::loadFeature(FILE * file){
  unsigned char type;
  fread(&type,sizeof(unsigned char),1,file);
  if(type==BASIC_FEATURE){
    Feature * feat=new Feature(type);
    if(feat->load(file))return feat;
    delete feat;
    return 0;
  }
  else if(type==IMAGE_FEATURE){
    Feature *feat=new ImageFeature();
    if(feat->load(file))return feat;
    delete feat;
    return 0;
  }
  else if(type==SIFT_FEATURE){
    SiftFeature * feat=new SiftFeature();
    if(feat->load(file))return feat;
    delete feat;
    return 0;
  }
  return 0;
}

std::vector<Feature*> getFeatures(Image8 & image){
  vector<Feature *> featuresi;
  //currently using both feature types
  int hw=4;

  //If we have an alpha channel, use it as a mask to filter out bad points.
  nacb::Image8 mask;
  if(image.nchannels==4){
    mask = image.getChannel(3).boxFilter(hw);
  }
  else {
    mask = nacb::Image8(image.w, image.h, 1);
    mask = 255;
  }
  
  vector<Vec2f> corns=harris(image,1000,1);
  

  for(int j=0;j<corns.size();j++){
    int x=(int)corns[j].x;
    int y=(int)corns[j].y;
    if(x>=hw && y>=hw && x+hw<image.w && y+hw<image.h && mask(x, y)>128){
      Imagef sim;
      sim=image.subimage(x-hw,y-hw,2*hw+1,2*hw+1);
      featuresi.push_back(new ImageFeature(corns[j],sim));
    }
  }
  //features.push_back(featuresi);
  if(1){
    Imagef imf1;
    imf1=image;
    if(imf1.nchannels>=3)
      imf1=(imf1.getChannel(0)+imf1.getChannel(1)+imf1.getChannel(2))*0.33333;
    //imf1=imf1.resize(imf1.w*2,imf1.h*2);
    vector<SiftFeature> f=lowe_sift(imf1);
    for(int j=0;j<f.size();j++){
      if(mask((int)f[j].point.x, (int)f[j].point.y)<128)continue;

      SiftFeature * sf=new SiftFeature();
      *sf=f[j];
      featuresi.push_back(sf);
    }
  }
  return featuresi;
}


vector<Matrix> getMatchedPoints(const vector<Feature *> & f1,
				const vector<Feature *> & f2,
				const vector<int> & matches){
  vector<Matrix> xs;

  assert(matches.size() == f1.size());

  int good=0;
  for(int k=0; k<matches.size(); k++){
    assert(matches[k]<=(int)f2.size());

    good+=(matches[k]>=0);
  }
  
  if(good<0)return xs;

  //check for homography as well.
  Matrix x1(2,good);
  Matrix x2(2,good);
  int ind=0;
  for(int k=0; k<matches.size(); k++){
    if(matches[k]>=0){
      x1(0,ind) = f1[k]->point.x;
      x1(1,ind) = f1[k]->point.y;
      x2(0,ind) = f2[matches[k]]->point.x;
      x2(1,ind) = f2[matches[k]]->point.y;
      ind++;
    }
  }
  xs.push_back(x1);
  xs.push_back(x2);
  return xs;
}
