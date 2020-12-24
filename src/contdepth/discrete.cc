#include "discrete.h"
#include "hermite.h"

#include "lbfgs.h"

#include "autorecon/stereo/conversions.h"
#include "autorecon/stereo/simplesequence.h"

#include "autorecon/recon_geometry.h"
#include "autorecon/recon_globals.h"

#include <nmath/vec3.h>
#include <nmath/vec2.h>

#include "MRF2.1/BP-S.h"


using namespace nacb;


const double kAlphaC = 10;
const double kEta = 0.02;
const float kMaxCost = 10.0;


float image_cost(const Vec3<int> & t1,
		 const Vec3<int> & t2){
  return kAlphaC/(kAlphaC + (t1 - t2).len());
}

double get_filtered_cost(std::vector<float> & costs, float maxCost){
  double c = 0;
  int last = std::min(costs.size()/2, costs.size() - 1);

  std::sort(costs.begin(), costs.end());

  for(int ti=0; ti<=last; ti++)
    c += costs[ti];
  
  if(isnan(c))
    return maxCost;

  return c;
}


double smooth_func_abs(double dist,
		       double slope,
		       double * g1 = 0,
		       double * g2 = 0){
  if(g1){
    if(dist == 0)
      *g1 = *g2 = 0;
    else if(dist >0 && dist < kEta){
      *g1 =  slope;
      *g2 = -slope;
    }
    else if (dist < 0 && dist > -kEta){
      *g1 = -slope;
      *g2 =  slope;
    }
  }
  return std::min(fabs(dist), kEta);
}


double smooth_func_SSD(double dist,
		       double slope,
		       double * g1 = 0,
		       double * g2 = 0){
  if(g1){
    *g1 =  2.0*dist*slope*10000.;
    *g2 = -2.0*dist*slope*10000.;
  }
  return dist*dist*10000.;
}


template <int num, int denom>
double smooth_func_hermite(double dist,
			   double slope,
			   double * g1 = 0,
			   double * g2 = 0){
  static float f[5] = {0, kEta, 0, kEta, 0};
  double index = dist/kEta;
  double scale = double(num)/double(denom);

  if(index <= -1){
    if(g1)
      *g1 = *g2 = 0;
    return kEta*scale;
  }
  if(index >= 1){
    if(g1)
      *g1 = *g2 = 0;
    return kEta*scale;
  }
  if(g1){
    *g1 = hermite_interp(f, index + 2, 5, true) * slope * (scale/kEta);
    *g2 = -(*g1);
  }
  return hermite_interp(f, index + 2, 5) * scale;
}


void test_hermite_smooth(){
  printf("x = [");
  for(double x=-2*kEta; x<=2*kEta; x+=kEta/50.0){
    double g1, g2;
    double v = smooth_func_hermite<2, 1>(x, 1.0, &g1, &g2);
    double vup = smooth_func_hermite<2, 1>(x + 1e-5, 1.0, &g1, &g2);
    printf("%f %f\n", g1, (vup - v)/1e-5); 
    //smooth_func_hermite(x, 1.0, );
  }
  printf("];\n");
}



double continuous_cost(DiscreteDisparityMapInt & dmap,
		       int width, 
		       int height,
		       int nd, const float * cost,
		       const float * h_weights,
		       const float * v_weights,
		       double (* smooth_func)(double, double, double *, double *),
		       double * disp ,
		       double * grad = 0){
  double E_data = 0;
  double E_smooth = 0;
  double slope = dmap.labelToDisparitySlope();

  if(grad)
    memset(grad, 0, sizeof(double)*width*height);

  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){
      int pixel = y*width+x;

      E_data += hermite_interp(cost + pixel*nd, disp[pixel], nd);
      
      if(grad)
	grad[y*width + x] += hermite_interp(cost + pixel*nd, disp[pixel], nd, true);

      if(x < width){
	double dist = (dmap.labelToDisparity(disp[pixel]) - 
		       dmap.labelToDisparity(disp[pixel + 1]));
	double sval, g1, g2;
	
	if(grad)
	  sval = smooth_func(dist, slope, &g1, &g2);
	else
	  sval = smooth_func(dist, slope, 0, 0);

	E_smooth += h_weights[pixel]*sval;

	if(grad){
	  double factor = h_weights[pixel];
	  grad[pixel] += factor*g1;
	  grad[pixel + 1] += factor*g2;
	}
      }
      if(y + 1 < height){
	double dist = (dmap.labelToDisparity(disp[pixel]) - 
		       dmap.labelToDisparity(disp[pixel + width]));
	double sval, g1, g2;
	
	if(grad)
	  sval = smooth_func(dist, slope, &g1, &g2);
	else
	  sval = smooth_func(dist, slope, 0, 0);

	E_smooth += v_weights[pixel]*sval;

	if(grad){
	  double factor = v_weights[pixel];
	  grad[pixel] += factor*g1;
	  grad[pixel + width] += factor*g2;
	}
      }
    }
  }
  return E_smooth + E_data;
}


struct cont_func_t {

  int width, height, nd;
  const float * cost; 
  const float * h_weights;
  const float * v_weights;
  DiscreteDisparityMapInt & disp;

  double (* smooth_func)(double, double, double *, double *);

  cont_func_t(int _width, 
	      int _height, int _nd,
	      const float * _cost, 
	      const float * _h, const float * _v,
	      DiscreteDisparityMapInt & _disp, 
	      double (* _smooth)(double, double, double *, double *)) :
    width(_width), height(_height), 
    nd(_nd), cost(_cost), 
    h_weights(_h), v_weights(_v), 
    disp(_disp), smooth_func(_smooth) {

  }
  
  double func(int n, double * x, double * g){    
    double val = continuous_cost(disp, width, height, nd, cost, h_weights, v_weights, smooth_func, x, g);
    printf("%lf\n", val);
    return val;
  }

  static double s_func(int n,double * x,double * g,const void * d){
    return ((cont_func_t*)d)->func(n, x, g);
  }
};


void brute_search_test(int width, int height, int nd,
		       const float * cost, 
		       const float * h_weights,
		       const float * v_weights,
		       DisparityMap<float> & dispImage, 
		       double (* smooth_func)(double, double, double *, double *)){
  double totalDist = 0;

  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){
      double dcur = dispImage(x, y);
      double cmin = 10;
      double dmin = dcur;
      //double ccur = cost[(y*width + x)*nd + int(dcur)];
      
      bool print = false;//  && ccur < 10;
      if(print)
	printf("x = [");
      
      for(double d = dcur - 1; d <= dcur + 1; d += 0.01){
	double v = hermite_interp(cost + (y*width + x)*nd, d, nd);
	if(print){
	  printf("%f ", v);
	}
	if(v < cmin){
	  cmin = v;
	  dmin = d;
	}
      }
      if(print)
	printf("]\n");
      
      dispImage(x, y) = dmin;
      double a = std::min(fabs(dmin - floor(dmin)), fabs(dmin - ceil(dmin)));
      totalDist += a;
    }    
  }
}


/**
   Yet another steepest descent optimization routine.
*/
void solve_continuous_descent(int width, int height, int nd,
			      const float * cost, 
			      const float * h_weights,
			      const float * v_weights,
			      DiscreteDisparityMapInt & disparityMapInt, 
			      DisparityMap<float> & dispImage, 
			      double (* smooth_func)(double, double, double *, double *)){
  const int maxIterations = 10;
  const bool checkGradient = false;
  double C = 0;
  double sc = 1.0;

  Matrix disp = disparityMapInt.getFlatMatrix();
  Matrix grad(width*height, 1);

  for(int t=0; t<maxIterations + checkGradient; t++){
    printf("Evaluating cost (it: %d)\n", t);
    C = continuous_cost(disparityMapInt, width, height, nd,
			cost, h_weights, v_weights, 
			smooth_func,
			disp.data, grad.data);
    printf("Cost is %lf\n", C);

    if(checkGradient){
      Matrix gfd(width*height, 1);
      for(int i=0; i<width*height; i++){
	if(fabs(grad[i])<1e-3)
	  continue;

	double cback = disp[i];
	disp[i] = disp[i] + 1e-4;
	double Cup = continuous_cost(disparityMapInt, width, height, nd,
				     cost, h_weights, v_weights, 
				     smooth_func,
				     disp.data, grad.data);
	printf("%f %f\n", grad[i], (Cup - C)/1e-4);
	disp[i] = cback;
      }
    }

    double s = sc;
    bool found = false;
    while(s >= 1e-16){
      Matrix dnew = disp.copy();
      dnew -= grad*s;
      
      double Cnew = continuous_cost(disparityMapInt, width, height, nd,
				    cost, h_weights, v_weights, 
				    smooth_func,
				    dnew.data);
      printf(" %lf s", Cnew);
      if(Cnew < C - 1e-4){
	sc = s*2.0;
	disp = dnew;
	found = true;
	break;
      }
      s /= 2;
    }
    if(!found)break;
  }
}


DisparityMap<float> solve_continuous(int width, int height, int nd,
				     const float * cost, 
				     const float * h_weights,
				     const float * v_weights,
				     DiscreteDisparityMapInt & disparityMapInt, 
				     double (* smooth_func)(double, double, double *, double *)){
  Matrix disp;
  bool bruteSearchTest = false;

  DisparityMap<float> dispImage(disparityMapInt);

  if(bruteSearchTest){
    brute_search_test(width, height, nd, cost, 
		      h_weights, v_weights, dispImage, smooth_func);
    dispImage.getNormalizedDisparityImage().write("/tmp/disp-cont-init.png");
  }

  disp = dispImage.getFlatMatrix();


  Matrix l(width * height, 1), u(width * height, 1);
  cont_func_t cfunc(width, height, nd, cost, 
		    h_weights, v_weights, disparityMapInt, smooth_func);
  int * nbd = new int[width * height];

  memset(nbd, 0, sizeof(int) * width * height);
  
  lbfgs(cont_func_t::s_func, width * height, disp.data, l.data, u.data, nbd, &cfunc, 1e10, 1e-4);

  dispImage = disp;
  dispImage.getNormalizedDisparityImage().write("/tmp/disp-cont.png");

  delete [] nbd;

  return dispImage;
}


/**
   init_disparity_best
   \brief Compute the disparity (best for each pixel, ignoring smoothness)
*/
void init_disparity_best(const nacb::Image8 & foreground, float * cost, DiscreteDisparityMapInt & disparityMap){
  int nd = disparityMap.nd;

  for(int y=0; y<foreground.height; y++){
    for(int x=0; x<foreground.width; x++){
      int best = 0;
      int pixel = y*foreground.width + x;

      if(!foreground(x, y)){
	disparityMap(x, y) = nd;
	continue;
      }
      
      for(int di=0; di<nd; di++){
	if(isnan(cost[pixel*nd + di]))
	  continue;
	
	if(cost[pixel*nd + di] < cost[pixel*nd + best]){
	  best = di;
	}	     
      }
      
      disparityMap(x, y) = best;
    }
  }
}



void cost_spatial_smoothness(const nacb::Image8 & image, 
			     const nacb::Image8 & foreground,
			     float * max_weights,
			     float * h_weights,
			     float * v_weights,
			     float smooth_weight,
			     float kSpatialEps = 0.1){

  int width = image.width, height = image.height;

  memset(max_weights, 0, sizeof(*max_weights)*width*height);

  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){
      nacb::Vec3<int> color = color_sample(image, x, y);
      
      if(x + 1 < width){
	Vec3<int> neigh = color_sample(image, x+1, y);
	float dist = 1.0/((color - neigh).len() + kSpatialEps);
	
	max_weights[y*width + x] += dist;
	max_weights[y*width + x + 1] += dist;
	h_weights[y*width + x] = dist;
      }
      if(y + 1 < height){
	Vec3<int> neigh = color_sample(image, x, y+1);
	float dist = 1.0/((color - neigh).len() + kSpatialEps);
	
	max_weights[y*width + x] += dist;
	max_weights[(y + 1)*width + x] += dist;
	v_weights[y*width + x] = dist;
      }
    }
  }
  
  // Normalize the smoothness cost
  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){
      int nneigh = (x > 0) + (x + 1 < width) + (y > 0) + (y + 1 < height);
      
      max_weights[y*width + x] /= nneigh;
      
      h_weights[y*width + x] /= max_weights[y*width + x];
      v_weights[y*width + x] /= max_weights[y*width + x];
      
      h_weights[y*width + x] *= smooth_weight;
      v_weights[y*width + x] *= smooth_weight;
      
      if(!foreground(x, y)){
	h_weights[y*width + x] = 0;
	v_weights[y*width + x] = 0;
	
	if(x > 0)
	  h_weights[y*width + (x-1)] = 0;
	if(y > 0)
	  v_weights[(y-1)*width + x] = 0;
      }
    }
  }
}


/**
   Parameters: kSigmaDepth, hwin, kMaxCost
 */
void cost_data(std::vector<SeqType> & seqs, int i, 
	       const std::vector<std::vector<EpipolarLine> > & epi_lines,
	       std::vector<DiscreteDisparityMapInt> & disparityMaps,
	       int width, int height, int nd, float * dmap, bool initialize, float * cost, int nclosest = -1){  
  const double kSigmaDepth = 3.0;

  printf("Getting disparity:\n");
  for(int di=0; di<nd; di++){
    if(di % 15 == 0) {
      printf(" %d ", di);
      fflush(stdout);
    }
      
    float d = dmap[di];

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	//int slice = di*width*height;
	  
	if((seqs[i].image.nchannels == 4) && (seqs[i].image(x, y, 3) < 20))
	  continue;
	  
	
	std::vector<float> costs;

	int ti_low, ti_high;
	get_closest_range(nclosest, seqs.size(), i, ti_low, ti_high);

	for(int ti=ti_low; ti<ti_high; ti++){
	  if(ti == i)continue;
	  
	  // Map the pixel to the neighboring images and compute difference.
	  Vec2<float> neighbor = epi_lines[i][ti](x, y, d);
	  float cost = image_cost(color_sample(seqs[i].image, x, y),
				  color_sample(seqs[ti].image, neighbor.x, neighbor.y));
	  float geomCost = 1.0;
	    
	  if(!initialize){
	    float neighborDisp;

	    if((neighborDisp = disparityMaps[ti].bilinear(neighbor.x, neighbor.y)) >= 0){
	      Vec2<float> back = epi_lines[ti][i](neighbor.x, neighbor.y, neighborDisp);
	      Vec2<float> diff = back - Vec2f(x, y);
	      double minDist = diff.dot(diff);

	      int hwin = 2;
	      int nx = (int)round(neighbor.x), ny = (int)round(neighbor.y);
		
	      for(int yy = std::max(0, ny - hwin); yy <= std::min(height - 1, ny + hwin); yy++){
		for(int xx = std::max(0, nx - hwin); xx <= std::min(width - 1, nx + hwin); xx++){
		  back = epi_lines[ti][i](neighbor.x, neighbor.y, disparityMaps[ti].getDisparity(xx, yy));
		  diff = back - Vec2f(x, y);
		  minDist = std::min(minDist, diff.dot(diff));
		}
	      }
	      geomCost = exp(-(minDist/(2.0*kSigmaDepth*kSigmaDepth)));
	    }
	    else
	      geomCost = 0.0f;
		
	  }
	  costs.push_back(cost * geomCost);
	}

	cost[(y*width + x)*nd + di] = get_filtered_cost(costs, kMaxCost);
      }
    }      
  }
  printf("\n");

  // Perform final normalization.
  for(int pi=0; pi<width*height; pi++){
    float maxv = 0;

    for(int di=0; di<nd; di++)
      maxv = std::max(maxv, cost[pi*nd + di]);
        
    for(int di=0; di<nd; di++)
      cost[pi*nd + di] = 1.0 - cost[pi*nd + di]/maxv;
  }
}


void solve_disparity(std::vector<SeqType> & seqs, 
		     float dmin, float dmax, int nd,
		     std::vector<DiscreteDisparityMapInt> & disparityMaps,
		     double smoothWeight,
		     int nclosest,
		     int downsample){
  std::vector<std::vector<EpipolarLine> > epi_lines(seqs.size());

  // Initialize the disparity maps for each image (e.g., sequence).
  int width = seqs[0].image.w;
  int height = seqs[0].image.h;

  float * cost = new float[width*height*nd];
  float * smooth_cost = new float[nd*nd];

  float * dmap = new float[nd];
  float * h_weights = new float[width*height];
  float * v_weights = new float[width*height];
  float * max_weights = new float[width*height];

  for(int di=0; di<nd; di++)
    dmap[di] = di/float(nd - 1)*(dmax - dmin) + dmin;
  

  for(int i=0; i<nd; i++){
    for(int j=0; j<nd; j++){
      smooth_cost[i*nd + j] = std::min(fabs(double(dmap[i] - dmap[j])), kEta);
    }
  }

  for(int i=0; i<(int)seqs.size(); i++){
    for(int j=0; j<(int)seqs.size(); j++){
      epi_lines[i].push_back(EpipolarLine(seqs[i].A, seqs[i].E,
					  seqs[j].A, seqs[j].E));
    }
  }

  bool initialize = (disparityMaps.size() == 0);

  if(initialize){
    for(int i=0; i<(int)seqs.size(); i++)
      disparityMaps.push_back(DiscreteDisparityMapInt(width, height, dmin, dmax, nd));
  }


  for(int i=0; i<(int)seqs.size(); i++){
    printf("on sequence %d\n", i);
    
    for(int pi=0; pi<nd*width*height; pi++)
      cost[pi] = kMaxCost;

    nacb::Image8 foreground;
    foreground = seqs[i].image.getChannel(3) > 128;
    
    cost_spatial_smoothness(seqs[i].image, foreground, max_weights, h_weights, v_weights, smoothWeight);
    cost_data(seqs, i, epi_lines, disparityMaps, width, height, nd, dmap, initialize, cost, nclosest);
        
    
    // Initialize the disparity map with best label.
    if(initialize)
      init_disparity_best(seqs[i].image, cost, disparityMaps[i]);

    DataCost dataCost(cost);
    SmoothnessCost smoothCost(smooth_cost, h_weights, v_weights);
    EnergyFunction energy(&dataCost, &smoothCost);
    BPS bps(width, height, nd, &energy);
    float tm;
    
    bps.initialize();
    bps.clearAnswer();
    
    for(int y=0; y<height; y++)
      for(int x=0; x<width; x++){
	bps.setLabel(y*width + x, (int)disparityMaps[i](x, y));
      }
      
    printf("optimizing..., %f = %f + %f\n", bps.totalEnergy(), bps.dataEnergy(), bps.smoothnessEnergy());
    
    for(int it=0; it<4; it++){
      bps.optimize(1, tm);
      printf("          ..., %f = %f + %f  (took:%f)\n", 
	     bps.totalEnergy(), bps.dataEnergy(), bps.smoothnessEnergy(), tm);
    }

    // Pull out the labels.
    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	int lab = bps.getLabel(y*width + x);
	
	disparityMaps[i](x, y) = lab;	  
      }
    }
    
    std::string snum = (boost::format("-%04d") % i).str();
    std::string dir = "/tmp/";
    nacb::Imagef depth = disparityMaps[i].getDepth(foreground, 0);

    {
      nacb::Imagef depthFull = depth.resize(depth.w << downsample, depth.h << downsample);
      depthFull.write((dir + "depth" + snum + ".png").c_str());
      writeRawDepth((dir + "depth" + snum + ".rfi").c_str(), depthFull);
    }
    disparityMaps[i].getNormalizedDisparityImage().write((dir + "disp" + snum + ".png").c_str());
    
    bool continuous = false;
    if(!initialize && continuous){
      printf("Doing continuous optimization.\n");
      DisparityMap<float> contMap = solve_continuous(width, height, nd, cost,
						     h_weights, v_weights,
						     disparityMaps[i], smooth_func_hermite<60, 1>);
      
      depth = contMap.getDepth(foreground, 0);
      depth.write((dir + "dcont" + snum + ".png").c_str());
      writeRawDepth((dir + "dcont" + snum + ".rfi").c_str(), depth);
    }
  }

  delete [] cost;
  delete [] h_weights;
  delete [] v_weights;
  delete [] dmap;
  delete [] smooth_cost;
}




/**
   return an offset vector
*/
std::vector<DiscreteFlowMap> solve_flow(std::vector<SeqType> & seqs, 
					std::vector<DiscreteDisparityMapInt> & disparityMaps,
					float dmax, int nx, int ny, int nz){
  std::vector<DiscreteFlowMap> flowMaps;
  int width = seqs[0].image.width;
  int height = seqs[0].image.height;

  for(int i=0; i<(int)seqs.size(); i++)
    flowMaps.push_back(DiscreteFlowMap(width, height, dmax, dmax, dmax, nx, ny, nz));

  int nlabels = flowMaps[0].getNumLabels();

  float * cost = new float[nlabels*width*height];
  float * smooth_cost = new float[nlabels*nlabels];
  float * h_weights = new float[width*height];
  float * v_weights = new float[width*height];
  float * max_weights = new float[width*height];


  // Get the smoothness cost, by unpacking the labels.
  for(int i=0; i<nlabels; i++){
    Vec3f di = flowMaps[0].getOffsetForLabel(i);
    
    for(int j=0; j<nlabels; j++){
      Vec3f dj = flowMaps[0].getOffsetForLabel(j);
      
      float dist = (di - dj).len();
      smooth_cost[i*nlabels + j] = dist; // FIXME: truncated?
      smooth_cost[j*nlabels + i] = dist;
    }
  }


  // We need to compute the cost volume and 
  for(int si=0; si<(int)seqs.size(); si++){
    nacb::Image8 foreground;
    foreground = seqs[si].image.getChannel(3) > 128;

    for(int i=0; i<nlabels*width*height; i++)
      cost[i] = kMaxCost;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	if((seqs[si].image.nchannels == 4) && (seqs[si].image(x, y, 3) < 20))
	  continue;

	float z = 1.0/disparityMaps[si].getDisparity(x, y);
	Vec3f wrld = backProject(seqs[si].P, float(x), float(y), z);

	for(int li = 0; li<nlabels; li++){
	  Vec3f offs = flowMaps[si].getOffsetForLabel(li);
	  std::vector<float> costs;
	      
	  for(int ti=0; ti<(int)seqs.size(); ti++){
	    if(ti == si)continue;
	    
	    float scaleOffset = (ti - si);
	    Vec3f worldOffset = wrld + offs*scaleOffset;
	    
	    Vec2f proj = project(seqs[ti].P, worldOffset);
	    float cost = image_cost(color_sample(seqs[si].image, x, y),
				    color_sample(seqs[ti].image, proj.x, proj.y));
	    float geomCost = 1.0;
#warning  Geometric cost is incomplete.
	    // if(!initialize) ...
	    costs.push_back(cost * geomCost);
	  }

	  cost[(y*width + x)*nlabels + li] = get_filtered_cost(costs, kMaxCost);	  
	}
      }
    }

    // Perform final normalization.
    for(int pi=0; pi<width*height; pi++){
      float maxv = 0;
      for(int di=0; di<nlabels; di++)
	maxv = std::max(maxv, cost[pi*nlabels + di]);
      
      for(int di=0; di<nlabels; di++)
	cost[pi*nlabels + di] = 1.0 - cost[pi*nlabels + di]/maxv;
    }

    cost_spatial_smoothness(seqs[si].image, foreground, max_weights, h_weights, v_weights, kSmoothWeight/5);

    DataCost dataCost(cost);
    SmoothnessCost smoothCost(smooth_cost, h_weights, v_weights);
    EnergyFunction energy(&dataCost, &smoothCost);
    BPS bps(width, height, nlabels, &energy);
    
    bps.initialize();
    bps.clearAnswer();
    

    // Initialize the labels to be zero.
    for(int y=0; y<height; y++)
      for(int x=0; x<width; x++)
	bps.setLabel(y*width + x, flowMaps[si].getZeroLabel());
    
    printf("optimizing..., %f = %f + %f\n", bps.totalEnergy(), bps.dataEnergy(), bps.smoothnessEnergy());
    
    for(int it=0; it<4; it++){
      float tm;
      bps.optimize(1, tm);

      std::cout << "Done optimizing:" << bps.totalEnergy() << "="
		<< bps.dataEnergy() << " + " << bps.smoothnessEnergy() << "\n";
    }

    nacb::Imagef dlen(width, height, 1);
    float maxlen = 1e-4;

    for(int y=0; y<height; y++){
      for(int x=0; x<width; x++){
	int label = bps.getLabel(y*width + x);
	flowMaps[si](x, y) = label;

	Vec3f dj = flowMaps[si].getOffset(x, y);
	dlen(x, y) = dj.len();

	maxlen = std::max(dlen(x, y), maxlen);
      }
    }
    dlen *= (1.0/maxlen);
    dlen.write("/tmp/dlen.png");
  }
  
  delete [] cost;
  
  return flowMaps;
}



