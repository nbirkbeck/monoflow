#ifndef CONTDEPTH_UTILS_H
#define CONTDEPTH_UTILS_H

#include <nimage/image.h>
#include "DisparitySequence.h"
#include "DisparityMap.h"
#include "DisparityOffset.h"
#include "FlowMap.h"
#include "EpipolarLine.h"

#include <boost/function.hpp>
#include <boost/bind.hpp>

typedef boost::function<nacb::Image8 (int i, int j)> warp_by_index_t;

struct image_stats_t {
  enum flags_t {
    NONE = 0,
    MEAN = 0x1,
    MIN = 0x2,
    MAX = 0x4,
    MEDIAN = 0x8,
    ALL = 0xFFFF
  }flags;

  double mn;
  double mx;
  double mean;
  double median;

  image_stats_t() {
    mean = 0;
    mn = 0;
    mx = 0;
    median = 0;
    flags = NONE;
  }
  
  image_stats_t(double _mn,
		double _mx,
		double _mean,
		double _median) : mn(_mn), mx(_mx), mean(_mean), median(_median) 
  {
    flags = ALL;
  }

  image_stats_t& operator+=(const image_stats_t& other) {
    mean += other.mean;
    mn = std::min(mn, other.mn);
    mx = std::max(mx, other.mx);
    median += other.median;
    return *this;
  }
  
  image_stats_t& operator/=(double d) {
    mean /= d;
    median /= d;
    return *this;
  }

  void print(std::ostream & o) const {
    if (flags == MEAN)
      o << mean;
    else
      o << "[" << mean << "," << median << "," << mn << "," << mx << "]";
  }
};


inline std::ostream& operator<<(std::ostream & o, const image_stats_t & stats){
  stats.print(o);
  return o;
}


image_stats_t image_compare(const nacb::Image8 & mask,
			    const nacb::Image8 & image1,
			    const nacb::Image8 & image2);

image_stats_t image_compare(const nacb::Imagef & mask,
			    const nacb::Imagef & image1,
			    const nacb::Imagef & image2);

void get_closest_range(const int nclosest, const int seqs, int i, int & lo, int & hi);


std::vector<Imagef> resize_images(const std::vector<Imagef> & images, int w, int h);

void save_depth_and_flow(const std::vector<Imagef> & flows,
			 const std::string & depthBase,
			 const std::string & flowBase,
			 bool UseDisparity = true);

inline nacb::Vec3<int> color_sample(const nacb::Image8 & image, float x, float y){
  unsigned char res[3];

  if(x <= 0 || y <= 0 || x>=image.w || y>=image.h){
    static int r = 0;
    r += 2048;
    return nacb::Vec3<int>(r, r, r);
  }

  image.bilinear(x, y, res, 3);
  return nacb::Vec3<int>(res[0], res[1], res[2]);
}

enum warp_alpha_t {
  WARP_ALPHA_NONE,
  WARP_ALPHA_SOURCE,
  WARP_ALPHA_DEST,
  WARP_ALPHA_WHITE,
  WARP_ALPHA_DEFAULT = WARP_ALPHA_SOURCE
};

// Separate disparity and flow.
template <class DisparityOffsetType>
nacb::Image8 warp_image(const nacb::Image8 & image,
			const DisparityOffsetType & warp,
			const nacb::Imagef & disp,
			const nacb::Imagef & flow,
			const warp_alpha_t & warp_alpha = WARP_ALPHA_DEFAULT);

nacb::Image8 warp_image(const nacb::Image8 & image,
			const EpipolarLine & eline,
			const nacb::Imagef & disp,
			const warp_alpha_t & warp_alpha = WARP_ALPHA_DEFAULT);

inline nacb::Image8 warp_image(const std::vector<SeqType> &seqs,
			       const std::vector<Imagef> & disps,
			       const std::vector<Imagef> & flows,
			       int i, int j,
			       const warp_alpha_t & warp_alpha = WARP_ALPHA_DEFAULT) {
  if (i == j) return seqs[i].image.copy();
  
  DisparityOffset warp(seqs[i].A, seqs[i].E,
		       seqs[j].A, seqs[j].E, double(j - i));

  nacb::Image8 image = warp_image(seqs[j].image, warp, disps[i], flows[i],
				  (warp_alpha == WARP_ALPHA_DEST) ? 
				  WARP_ALPHA_SOURCE : warp_alpha);

  if (warp_alpha == WARP_ALPHA_DEST && image.nchannels == 4) 
    image.setChannel(3, seqs[i].image.getChannel(3));

  return image;
}


inline nacb::Image8 warp_image(const std::vector<SeqType> &seqs, 
			       const std::vector<nacb::Imagef> & disps,
			       int i, int j,
			       const warp_alpha_t & warp_alpha = WARP_ALPHA_DEFAULT) {
  EpipolarLine eline(seqs[i].A, seqs[i].E,
		     seqs[j].A, seqs[j].E);
  
  nacb::Image8 image = warp_image(seqs[j].image, eline, disps[i],
				  (warp_alpha == WARP_ALPHA_DEST) ? 
				  WARP_ALPHA_SOURCE : warp_alpha);

  if (warp_alpha == WARP_ALPHA_DEST && image.nchannels == 4) 
    image.setChannel(3, seqs[i].image.getChannel(3));

  return image;
}


template <class DisparityOffsetType>
nacb::Image8 warp_image(int w, int h, 
			const std::vector<nacb::Imagef> & images,
			const std::vector<DisparityOffsetType> & warps,
			const nacb::Imagef & data,
			int j);

nacb::Image8 warp_image(const std::vector<SeqType> &seqs, 
			const std::vector<DiscreteDisparityMapInt> & disps,
			int i, int j);

nacb::Image8 warp_image(const std::vector<SeqType> &seqs, 
			const std::vector<DiscreteDisparityMapInt> & disps,
			const std::vector<DiscreteFlowMap> & flows,
			int i, int j);

void write_all_warped_images(const warp_by_index_t & warp_func,
			     int nseqs,
 			     const char * prefix,
			     int nclosest = -1);

inline void write_all_warped_images(const std::vector<SeqType> &seqs, 
				    const std::vector<DiscreteDisparityMapInt> & disps,
				    const char * prefix,
				    int nclosest = -1){
  // Cast overloaded function
  nacb::Image8 (* warp_image_func)(const std::vector<SeqType> &seqs, 
				   const std::vector<DiscreteDisparityMapInt> & disps,
				   int i, int j) = warp_image;

  warp_by_index_t warper = boost::bind(warp_image_func, seqs, disps, _1, _2);

  write_all_warped_images(warper, seqs.size(), prefix, nclosest);
}


inline void write_all_warped_images(const std::vector<SeqType> &seqs, 
				    const std::vector<Imagef> & disps,
				    const std::vector<Imagef> & flows,
				    const char * prefix,
				    int nclosest = -1,
				    const warp_alpha_t & warp_alpha = WARP_ALPHA_DEFAULT){
  // Cast overloaded function
  nacb::Image8 (* warp_image_func)(const std::vector<SeqType> & seqs, 
				   const std::vector<Imagef> & disps,
				   const std::vector<Imagef> & flows,
				   int i, int j, const warp_alpha_t &) = warp_image;
  
  warp_by_index_t warper = boost::bind(warp_image_func, seqs, disps, flows, _1, _2, warp_alpha);

  write_all_warped_images(warper, seqs.size(), prefix, nclosest);
}


inline void write_all_warped_images(const std::vector<SeqType> & seqs, 
				    const std::vector<Imagef> & disps,
				    const char * prefix,
				    int nclosest = -1,
				    const warp_alpha_t & warp_alpha = WARP_ALPHA_DEFAULT){
  // Cast overloaded function
  nacb::Image8 (* warp_image_func)(const std::vector<SeqType> & seqs, 
				   const std::vector<Imagef> & disps,
				   int i, int j, const warp_alpha_t &) = warp_image;
  
  warp_by_index_t warper = boost::bind(warp_image_func, seqs, disps, _1, _2, warp_alpha);

  write_all_warped_images(warper, seqs.size(), prefix, nclosest);
}



void write_all_warped_images(std::vector<SeqType> &seqs, 
			     const std::vector<DiscreteDisparityMapInt> & disps,
			     const std::vector<DiscreteFlowMap> & flows,
			     const char * prefix);

nacb::Imagef disparity_to_depth(const nacb::Imagef & disp);

inline nacb::Imagef depth_to_disparity(const nacb::Imagef & disp){
  return disparity_to_depth(disp);
}

image_stats_t flow_compare(const nacb::Image8 & image,
			   const nacb::Imagef & truth,
			   const nacb::Imagef & flow);

image_stats_t  disp_compare(const nacb::Image8 & image,
			    const nacb::Imagef & truth, 
			    const nacb::Imagef & disp);

nacb::Imagef weight_from_alpha(const nacb::Image8 & image, int blur = 3);


// 3D flow conversions
nacb::Imagef project_flow(const nacb::Matrix & A, const nacb::Matrix & E,
			  const nacb::Imagef & dispFlow);


#endif //CONTDEPTH_UTILS_H
