#ifndef DISCRETE_H
#define DISCRETE_H

#include "utils.h"
#include "DisparityMap.h"
#include "DisparitySequence.h"
#include "FlowMap.h"
#include "EpipolarLine.h"

#include <vector>
#include <nmath/vec3.h>

static const double kSmoothWeight = 250.0;//300.0;


inline nacb::Vec3<int> color_sample(const nacb::Image8 & image, int x, int y){
  return nacb::Vec3<int>(image(x, y, 0), image(x, y, 1), image(x, y, 2));
}

void solve_disparity(std::vector<SeqType> & seqs, 
		     float dmin, float dmax, int nd,
		     std::vector<DiscreteDisparityMapInt> & disparityMaps,
		     double smoothWeight,
		     int nclosest,
		     int downsample);

std::vector<DiscreteFlowMap> solve_flow(std::vector<SeqType> & seqs, 
					std::vector<DiscreteDisparityMapInt> & disparityMaps,
					float dmax, int nx, int ny, int nz);
  

#endif // DISCRETE_H
