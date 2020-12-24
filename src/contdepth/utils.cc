#include "utils.h"
#include "EpipolarLine.h"

#include <nmath/vec3.h>
#include <nmath/vec2.h>

using namespace nacb;


image_stats_t image_compare(const nacb::Imagef & mask,
			    const nacb::Imagef & image1,
			    const nacb::Imagef & image2){
  assert(mask.nchannels == 1);

  int nch = std::min(3, image1.nchannels);

  double diff = 0;
  double ndiff = 0;
  for(int y=0; y<image1.h; y++){
    for(int x=0; x<image1.w; x++){

      if(mask(x, y) > 0){
	for(int k=0; k<nch; k++){
	  double a = image1(x, y, k);
	  double b = image2(x, y, k);
	  a -= b;
	  diff += a*a * mask(x, y); 
	}	
	ndiff += mask(x, y)*nch;
      }
    }
  }
  image_stats_t stats;
  stats.mean = sqrt(diff/ndiff);
  stats.flags = image_stats_t::MEAN;
  return stats;
}


image_stats_t image_compare(const nacb::Image8 & mask,
			    const nacb::Image8 & image1,
			    const nacb::Image8 & image2){
  assert(mask.nchannels == 1);

  int nch = std::min(3, image1.nchannels);

  double diff = 0;
  int ndiff = 0;
  for(int y=0; y<image1.h; y++){
    for(int x=0; x<image1.w; x++){

      if(mask(x, y)){
	for(int k=0; k<nch; k++){
	  double a = image1(x, y, k);
	  double b = image2(x, y, k);
	  a -= b;
	  diff += a*a; 
	}	
	ndiff++;
      }
    }
  }
  image_stats_t stats;
  stats.mean = sqrt(diff/ndiff);
  stats.flags = image_stats_t::MEAN;
  return stats;
}


void get_closest_range(const int nclosest, const int nseqs, int i, int & lo, int & hi){
  if (nclosest < 0) {
    lo = 0;
    hi = nseqs;
    return;
  }

  lo = std::max(0, i - nclosest);
  hi = std::min(nseqs, i + nclosest + 1);
  
  if(lo == 0)
    hi = std::min(2*nclosest + 1, nseqs);
  
  else if(hi == nseqs)
    lo = std::max(nseqs - 2*nclosest - 1, 0);
}


std::vector<Imagef> resize_images(const std::vector<Imagef> & images, 
				  int w, int h) {
  std::vector<Imagef> newImages;

  for(int i=0; i<(int)images.size(); i++)
    newImages.push_back(images[i].resize(w, h));

  return newImages;
}


void save_depth_and_flow(const std::vector<Imagef> & flows,
			 const std::string & depthBase,
			 const std::string & flowBase,
			 bool UseDisparity) {
  for(int i=0; i<(int)flows.size(); i++){
    if(depthBase.size()){
      std::string depthName = (boost::format(depthBase) % i).str();

      if (UseDisparity)
	disparity_to_depth(flows[i].getChannel(0)).write(depthName.c_str());
      else {
	flows[i].getChannel(0).write(depthName.c_str());
      }
    }

    if(flowBase.size()){
      std::string flowName = (boost::format(flowBase) % i).str();
      nacb::Imagef flow(flows[i].w, flows[i].h, 3);

      for(int k=0; k<3; k++)
	flow.setChannel(k, flows[i].getChannel(k + 1));
      
      flow.save(flowName.c_str());
    }
  }
}


nacb::Image8 warp_image(const nacb::Image8 & image,
			const EpipolarLine & eline,
			const nacb::Imagef & disp,
			const warp_alpha_t & warp_alpha) {
  int nchannels = warp_alpha == WARP_ALPHA_WHITE ? 4 : 
    ((warp_alpha == WARP_ALPHA_NONE) ? 3 : image.nchannels);

  nacb::Image8 warped(disp.w, disp.h, nchannels);

  for(int y=0; y<image.h; y++){
    for(int x=0; x<image.w; x++){
      nacb::Vec2f coord = eline(x, y, disp(x, y));
      nacb::Vec3<int> color = color_sample(image, coord.x, coord.y);
      
      for(int k=0; k<3; k++)
	warped(x, y, k) = color.data[k];

      if (nchannels == 4) {
	if (warp_alpha == WARP_ALPHA_WHITE)
	  warped(x, y, 3) = 255;
	else
	  warped(x, y, 3) = image.bilinear(coord.x, coord.y, 3);
      }
    }
  }
  return warped;
}


template <class DisparityOffsetType>
nacb::Image8 warp_image(const nacb::Image8 & image,
			const DisparityOffsetType & warp,
			const nacb::Imagef & disp,
			const nacb::Imagef & flow,
			const warp_alpha_t & warp_alpha) {
  int nchannels = warp_alpha == WARP_ALPHA_WHITE ? 4 : 
    ((warp_alpha == WARP_ALPHA_NONE) ? 3 : image.nchannels);
			
  nacb::Image8 result(disp.w, disp.h, nchannels);
 
  for(int y=0; y<image.h; y++){
    for(int x=0; x<image.w; x++){
      Vec3f offs(flow(x, y, 0), flow(x, y, 1), flow(x, y, 2));
      Vec2f warped = warp(x, y, disp(x, y), offs);
      Vec3<int> color = color_sample(image, warped.x, warped.y);

      for(int k=0; k<3; k++)
	result(x, y, k) = color[k];

      if (nchannels == 4) {
	if (warp_alpha == WARP_ALPHA_WHITE)
	  result(x, y, 3) = 255;
	else
	  result(x, y, 3) = image.bilinear(warped.x, warped.y, 3);
      }
    }	  
  }
  return result;
}


template <class DisparityOffsetType>
nacb::Image8 warp_image(int w, int h, 
			const std::vector<Imagef> & images,
			const std::vector<DisparityOffsetType> & warps,
			const nacb::Imagef & data,
			int j){ 
  nacb::Image8 image(w, h, 4);
 
  for(int y=0; y<h; y++){
    for(int x=0; x<w; x++){
      Vec3f offs(data(x, y, 1), data(x, y, 2), data(x, y, 3));
      Vec2f warped = warps[j](x, y, data(x, y, 0), offs);
      float color[4] = {0, 0, 0, 1};
      images[j].bilinear(warped.x, warped.y, color);
      
      for(int k=0; k<4; k++)
	image(x, y, k) = (unsigned char)std::max(0.f, std::min(255.f, color[k]*255));
    }
  }
  return image;
}


nacb::Image8 warp_image(const std::vector<SeqType> &seqs, 
			const std::vector<DiscreteDisparityMapInt> & disps,
			int i, int j){
  if(i == j)
    return seqs[i].image.copy();

  EpipolarLine epi(seqs[i].A, seqs[i].E, seqs[j].A, seqs[j].E);
  nacb::Image8 image(seqs[i].image.w, seqs[i].image.h, 4);
  
  for(int y=0; y<seqs[i].image.h; y++){
    for(int x=0; x<seqs[i].image.w; x++){
      Vec2f warped = epi(x, y, disps[i].getDisparity(x, y));
      Vec3<int> color = color_sample(seqs[j].image, warped.x, warped.y);
      
      for(int k=0; k<3; k++)
	image(x, y, k) = color.data[k];

      image(x, y, 3) = seqs[i].image.bilinear(x, y, 3);
    }
  }
  return image;
}


nacb::Image8 warp_image(const std::vector<SeqType> &seqs, 
			const std::vector<DiscreteDisparityMapInt> & disps,
			const std::vector<DiscreteFlowMap> & flows,
			int i, int j){
  if(i == j)
    return seqs[i].image.copy();

  EpipolarLine epi(seqs[i].A, seqs[i].E, seqs[j].A, seqs[j].E);
  nacb::Image8 image(seqs[i].image.w, seqs[i].image.h, 3);
  
  for(int y=0; y<seqs[i].image.h; y++){
    for(int x=0; x<seqs[i].image.w; x++){
      Vec3f point = backProject(seqs[i].P, float(x), float(y), 
				1.0f/disps[i].getDisparity(x, y)) + flows[i].getOffset(x, y);

      Vec2f p = project(seqs[j].P, point);
      Vec3<int> color = color_sample(seqs[j].image, p.x, p.y);
      
      for(int k=0; k<3; k++)
	image(x, y, k) = color.data[k];
    }
  }
  return image;
}


void write_all_warped_images(const warp_by_index_t & warp_func,
			     int nseqs,
 			     const char * prefix,
			     int nclosest){
  printf("writing all warped images\n");
  
  for(int i=0; i<nseqs; i++){
    printf("on seqeuence: %d\n", i);

    int low, high;
    get_closest_range(nclosest, nseqs, i, low, high);

    printf("range: %d %d\n",low, high);

    for(int j=low; j<high; j++){
      char fname[1024];
      snprintf(fname, sizeof(fname), prefix, i, j);

      warp_func(i, j).write(fname);
    }
  }
}


void write_all_warped_images(std::vector<SeqType> &seqs, 
			     const std::vector<DiscreteDisparityMapInt> & disps,
			     const std::vector<DiscreteFlowMap> & flows,
			     const char * prefix,
			     int nclosest){

  nacb::Image8 (* warp_image_func)(const std::vector<SeqType> & seqs, 
				   const std::vector<DiscreteDisparityMapInt> & disps,
				   const std::vector<DiscreteFlowMap> & flows,
				   int i, int j) = warp_image;

  write_all_warped_images(boost::bind(warp_image_func, seqs, disps, flows, _1, _2), 
			  seqs.size(), prefix, nclosest);
}


nacb::Imagef disparity_to_depth(const nacb::Imagef & disp){
  nacb::Imagef depth(disp.w, disp.h, 1);

  for(int y=0; y<disp.h; y++){
    for(int x=0; x<disp.w; x++){
      depth(x, y) = 1.0/disp(x, y);
    }
  }
  
  return depth;
}


image_stats_t flow_compare(const nacb::Image8 & image,
			   const nacb::Imagef & truth, const nacb::Imagef & flow){
  nacb::Imagef diff = truth - flow;
  std::vector<double> fdiffs;
  float minDiff = 1e8;
  float maxDiff = -1e8;
  
  double averageDistance = 0;

  for(int y=0; y<truth.h; y++){
    for(int x=0; x<truth.w; x++){
      if(image(x, y, 3)>100){
	Vec3f a(diff(x, y, 0), diff(x, y, 1), diff(x, y, 2));

	float d = a.len();
	fdiffs.push_back(fabs(d));	
	
	minDiff = std::min(minDiff, d);
	maxDiff = std::max(maxDiff, d);
	
	averageDistance += d;
      }
    }
  }
  averageDistance /= fdiffs.size();
  std::sort(fdiffs.begin(), fdiffs.end());
  return image_stats_t(minDiff, maxDiff, averageDistance, fdiffs[fdiffs.size()/2]);
}


image_stats_t disp_compare(const nacb::Image8 & image,
			   const nacb::Imagef & truth, const nacb::Imagef & disp){
  nacb::Imagef diff = truth - disp;
  
  std::vector<double> diffs;
  std::vector<double> fdiffs;
  float minDisp = 1e8, minTruth = 1e8, minDiff = 1e8;
  float maxDisp = -1e8, maxTruth = -1e8, maxDiff = -1e8;
  
  double averageDistance = 0;

  for(int y=0; y<truth.h; y++){
    for(int x=0; x<truth.w; x++){
      if(image(x, y, 3)>100){
	float d = diff(x, y);
	diffs.push_back(d);
	fdiffs.push_back(fabs(d));
	
	minTruth = std::min(minTruth, truth(x, y));
	maxTruth = std::max(maxTruth, truth(x, y));
	
	minDisp = std::min(minDisp, disp(x, y));
	maxDisp = std::max(maxDisp, disp(x, y));
	
	minDiff = std::min(minDiff, d);
	maxDiff = std::max(maxDiff, d);
	
	averageDistance += fabs(d);
      }
    }
  }
  averageDistance /= fdiffs.size();
  std::sort(fdiffs.begin(), fdiffs.end());

  //printf("truth: %f - %f,  disp: %f - %f\n", minTruth, maxTruth, minDisp, maxDisp);
  //printf("median distance: %lf\n", fdiffs[fdiffs.size()/2]);
  //printf("average: %f\n", averageDistance);
  
  return image_stats_t(minDiff, maxDiff, averageDistance, fdiffs[fdiffs.size()/2]);
}


nacb::Imagef weight_from_alpha(const nacb::Image8 & image, int blur){
  nacb::Imagef weights(image.getChannel(3)>1);
  weights *= 255;
  return weights.boxFilter(blur);
}



// 3D flow conversions
nacb::Imagef project_flow(const nacb::Matrix & A, const nacb::Matrix & E, 
			  const nacb::Imagef & dispFlow){
  //assert(disp.w == flow.w && disp.h == flow.h);
  assert(dispFlow.nchannels == 4);

  nacb::Imagef uv(dispFlow.w, dispFlow.h, 2);
  
  DisparityOffset doffs(A, E, A, E);

  for(int y=0; y<dispFlow.h; y++){
    for(int x=0; x<dispFlow.w; x++){
      Vec2d co = doffs(x, y, dispFlow(x, y ,0), 
		       Vec3d(dispFlow(x, y, 1),
			     dispFlow(x, y, 2),
			     dispFlow(x, y, 3)));
      uv(x, y, 0) = co.x - x;
      uv(x, y, 1) = co.y - y;
    }
  }
  return uv;
}


// Explicit instantiations for disparity or depth-based implementation.

template 
nacb::Image8 warp_image(const nacb::Image8 & ,
			const DisparityOffsetBase<true> &,
			const nacb::Imagef & ,
			const nacb::Imagef & ,
			const warp_alpha_t & warp_alpha);

template
nacb::Image8 warp_image(const nacb::Image8 & ,
			const DisparityOffsetBase<false> &,
			const nacb::Imagef & ,
			const nacb::Imagef & ,
			const warp_alpha_t & warp_alpha);

template
nacb::Image8 warp_image(int w, int h, 
			const std::vector<nacb::Imagef> & images,
			const std::vector<DisparityOffsetBase<true> > & warps,
			const nacb::Imagef & data,
			int j);

template
nacb::Image8 warp_image(int w, int h, 
			const std::vector<nacb::Imagef> & images,
			const std::vector<DisparityOffsetBase<false> > & warps,
			const nacb::Imagef & data,
			int j);
