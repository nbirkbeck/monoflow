#include <xmmintrin.h>
#include <nmath/matrix.h>
#include <vector>
#include <assert.h>

#include "flow.h"
#include "../recon_geometry.h"

using std::vector;
using nacb::Imagef;
using nacb::Matrix;

const static float flo_tag = 202021.25;

bool flow_save_crossfade(const char * basename,
			  const nacb::Imagef & im1,
			  const nacb::Imagef & im2,
			  const nacb::Imagef & d1,
			  const nacb::Imagef & d2,
			  int nt){
  std::vector<nacb::Imagef> results = flow_crossfade(im1, im2, d1, d2, nt);

  for(int i=0; i<(int)results.size(); i++){
    char name[2048];
    snprintf(name, 2048, basename, i); 
    results[i].write(name);
  }
  return true;
}

std::vector<nacb::Imagef> flow_crossfade(const nacb::Imagef & im1,
					 const nacb::Imagef & im2,
					 const nacb::Imagef & d1,
					 const nacb::Imagef & d2,
					 int nt){
  assert(im1.w == d1.w && im1.h == d1.h);
  assert(im2.w == d2.w && im2.h == d2.h);
  
  assert(d1.nchannels == 2);
  assert(d2.nchannels == 2);
  
  nt = std::max(nt, 3);

  std::vector<nacb::Imagef> results;

  for(int i=0; i<nt; i++){
    double t = double(i)/double(nt-1.0);

    //Want to warp back im2 a lot, e.g., (1.0-t)
    nacb::Imagef b = flow_warp_back(d1*(1.0-t), im2);

    //Want to only warp im1 forward a little (e.g., t)
    nacb::Imagef a = flow_warp_back(d2*t, im1);

    //Want to fade from im1_w to im2_w, so give more weight, (1.0-t) to a
    results.push_back(b*t + a*(1.0-t));
  }
  return results;
}

nacb::Matrix flow_fundamental_matrix(const nacb::Imagef & disp, const nacb::Image8 & mask,
				     Matrix * p1out, Matrix * p2out){
  int npts = 0;
    
  assert(disp.w == mask.w && disp.h == mask.h);
  assert(disp.nchannels == 2);

  for(int y=0; y<disp.h; y++)
    for(int x=0; x<disp.w; x++)
      if(mask(x, y))npts++;

  Matrix p1(2, npts);
  Matrix p2(2, npts);
  int ind = 0;

  for(int y=0; y<disp.h; y++){
    for(int x=0; x<disp.w; x++){
      if(mask(x, y)){
	p1(0, ind) = x;
	p1(1, ind) = y;
	
	p2(0, ind) = double(x) + disp(x, y, 0);
	p2(1, ind) = double(y) + disp(x, y, 1);
	ind++;
      }
    }
  }
  if(p1out)*p1out = p1;
  if(p2out)*p2out = p2;

  Matrix F = getFundamentalMatrix(p1, p2);
  return F;
}

bool flow_save_length(const nacb::Imagef & disp, const char * name){
  nacb::Imagef len;
  if(disp.nchannels==2)
    len = flow_length(disp);
  else
    len = disp;
  
  double n, x;
  len.getRange(n, x);
  
  return ((len - n)*(1.0/(x-n))).write(name);
}


nacb::Imagef flow_length(const nacb::Imagef & disp){
  assert(disp.nchannels == 2);

  nacb::Imagef len(disp.w, disp.h, 1);

  for(int y=0; y<disp.h; y++){
    for(int x=0; x<disp.w; x++){
      len(x, y) = sqrt(disp(x, y, 0)*disp(x, y, 0) + disp(x, y, 1)*disp(x, y, 1));
    }
  }
  return len;
}

//Middlebury file format.
bool flow_save_flo(const nacb::Imagef & image, const char * fname){
  if(image.nchannels != 2)return false;

  FILE * file = fopen(fname, "wb");
  if(!file)return false;

  float tag_true =  flo_tag;
  int dims[2] = {image.w, image.h};

  if(!fwrite(&tag_true, sizeof(float), 1, file)){
    fclose(file);
    return false;
  }

  if(fwrite(dims, sizeof(int), 2, file) != 2){
    fclose(file);
    return false;
  }
  
  fwrite(image.data, sizeof(float), dims[0]*dims[1]*2, file);
  fclose(file);
  return true;
}

//Middlebury file format
nacb::Imagef flow_read_flo(const char * fname){
  nacb::Imagef empty;
  FILE * file = fopen(fname, "rb");
  if(!file)return empty;
  float tag = 0;
  int dims[2];

  if(!fread(&tag, sizeof(float), 1, file) || tag != flo_tag){
    fclose(file);
    return empty;
  }

  if(fread(dims, sizeof(int), 2, file) != 2){
    fclose(file);
    return empty;
  }
  nacb::Imagef flo(dims[0], dims[1], 2);
  size_t r = fread(flo.data, sizeof(float), dims[0]*dims[1]*2, file);
  if ((int)r != dims[0]*dims[1]*2) {
    fprintf(stderr, "%s:%d expected %ld but read %d\n",
            __FILE__, __LINE__, r, dims[0]*dims[1]*2);
  }
  fclose(file);
  return flo;
}

/**
   Returns the compositional combination of two displacement fields U and update (upd)

   Notice that the compositional methods updates the flow by transforming the update 
   by looking up the displacement from existing displacement using newly displaced
   coordinates and then adding to this the updated displacement.
*/
nacb::Imagef flow_compose(const nacb::Imagef & U, const nacb::Imagef & upd){
  assert(U.w == upd.w && U.h == upd.h);
  assert(U.nchannels == 2 && upd.nchannels == 2);
  
  nacb::Imagef U2(U.w, U.h, 2);

  for(int y=0; y<U.h; y++){
    for(int x=0; x<U.w; x++){
      float xf = x + upd(x,y,0);
      float yf = y + upd(x,y,1);
      
      float dx = U.bilinear(xf, yf, 0);
      float dy = U.bilinear(xf, yf, 1);
      
      U2(x,y,0) = dx + upd(x,y,0);
      U2(x,y,1) = dy + upd(x,y,1);
    }
  }
  return U2;
}


nacb::Imagef flow_warp_back(const nacb::Imagef & disp, const nacb::Imagef & im1){
  nacb::Imagef warped(im1.w, im1.h, im1.nchannels);

  assert(disp.nchannels == 2);

  for(int y=0; y<im1.h; y++){    
    for(int x=0; x<im1.w; x++){
      float px = x + disp(x,y,0);
      float py = y + disp(x,y,1);
      for(int k=0; k<im1.nchannels; k++){
	warped(x,y,k) = im1.bilinear(px, py, k);
      }
    }
  }
  return warped;
}

nacb::Image8  flow_ensure_consistency(const nacb::Imagef & forward, const nacb::Imagef & backward, 
				      float forbackthresh, int badvalue){
  nacb::Image8 mask(forward.w, forward.h, 1);

  mask = 255;

  for(int y=0; y<forward.h; y++){
    for(int x=0; x<forward.w; x++){
      float x2 = std::max(0.f, std::min(float(x)+forward(x,y,0), (float)(forward.w-1)));
      float y2 = std::max(0.f, std::min(float(y)+forward(x,y,1), (float)(forward.h-1)));
      
      Vec2f boffs;
      backward.bilinear(x2, y2, boffs.data);
      float xback = x2 + boffs.x;
      float yback = y2 + boffs.y;
      
      if(fabs(xback-x)>=forbackthresh || fabs(yback-y)>=forbackthresh){	
	forward(x,y,0) = badvalue;
	forward(x,y,1) = badvalue;

	mask(x, y) = 0;
      }
    }
  }
  return mask;
}


//This provides much better results than flow_ensure_consistency...which is really additive.
nacb::Image8  flow_ensure_consistency_compose(const nacb::Imagef & forward, 
					      const nacb::Imagef & backward, 
					      float forbackthresh, int badvalue){
  nacb::Image8 mask(forward.w, forward.h, 1);

  mask = 255;

  nacb::Imagef composed = flow_compose(forward, backward);

  for(int y=0; y<forward.h; y++){
    for(int x=0; x<forward.w; x++){
      float x2 = composed(x, y, 0);
      float y2 = composed(x, y, 1);
            
      if(fabs(x2)>=forbackthresh || fabs(y2)>=forbackthresh){	
	forward(x,y,0) = badvalue;
	forward(x,y,1) = badvalue;

	mask(x, y) = 0;
      }
    }
  }
  return mask;
}

Imagef getSSD_sse(const Imagef & im1, const Imagef & im2,
		  int hwin,
		  int dx, int dy){
  int xstart = std::max(0, -dx);
  int xend   = std::min(im1.w, im2.w-dx);

  int ystart = std::max(0, -dy);
  int yend   = std::min(im1.h, im2.h-dy);

  int winsize = 2*hwin+1;

  Imagef result(im1.w, im1.h, 1);
  int nchannels = std::min(3, im1.nchannels);

  Imagef hpass(im1.w, im1.h, 1);

  hpass = hwin*nchannels+1;
  result = hwin*nchannels+1;


  //Horizontal pass
  for(int y=ystart; y<yend; y++){
    float sum = 0;
    for(int x=xstart; x<xstart+winsize; x++){
      for(int k=0; k<nchannels; k++){
	double diff = im1(x,y,k) - im2(x+dx,y+dy,k);
	sum += diff*diff;
      }
    }
    int x = 0;
    for(x=xstart+hwin; x<xend-hwin-1; x++){
      hpass(x,y) = sum;

      __m128 l1=_mm_loadu_ps(&im1(x-hwin,y,0));
      __m128 l2=_mm_loadu_ps(&im2(x-hwin+dx,y+dy,0));

      __m128 u1=_mm_loadu_ps(&im1(x+hwin+1,y,0));
      __m128 u2=_mm_loadu_ps(&im2(x+dx+hwin+1,y+dy,0));
      
      __m128 s1 = _mm_sub_ps(l1, l2);
      __m128 s2 = _mm_sub_ps(u1, u2);
      
      l1 = _mm_mul_ps(s1, s1);
      l2 = _mm_mul_ps(s2, s2);
      
      s1 = _mm_sub_ps(l2, l1);
      float unpacked[4];
      _mm_storeu_ps(unpacked, s1);
      sum += unpacked[0]+unpacked[1]+unpacked[2];
      hpass(x,y) = sum;
    }
    hpass(x,y) = sum;
  }

  float zero_f = 0.0f;
  
  //Vertical pass.
  for(int x=xstart; x<xend; x+=4){
    __m128 sum = _mm_load1_ps(&zero_f);
    for(int y=ystart; y<ystart+winsize; y++){
      __m128 d = _mm_loadu_ps(&hpass(x,y));
      sum = _mm_add_ps(sum, d);
    }
    
    int y = 0;
    for(y=ystart+hwin; y<yend-hwin-1; y++){
      _mm_storeu_ps(&result(x,y), sum);
      
      __m128 dl = _mm_loadu_ps(&hpass(x,y-hwin));
      __m128 du = _mm_loadu_ps(&hpass(x,y+hwin+1));
      
      sum = _mm_add_ps(_mm_sub_ps(du, dl), sum);
    }
    _mm_storeu_ps(&result(x,y), sum);
  }
  return result;  
}


Imagef getSSD(const Imagef & im1, const Imagef & im2,
	      int hwin,
	      int dx, int dy){
  int xstart = std::max(0, -dx);
  int xend   = std::min(im1.w, im2.w-dx);

  int ystart = std::max(0, -dy);
  int yend   = std::min(im1.h, im2.h-dy);

  int winsize = 2*hwin+1;

  Imagef hpass(im1.w, im1.h, 1);
  Imagef result(im1.w, im1.h, 1);
  int nchannels = std::min(3, im1.nchannels);

  hpass = hwin*nchannels+1;
  result = hwin*nchannels+1;
  //Horizontal pass
  for(int y=ystart; y<yend; y++){
    float sum = 0;
    for(int x=xstart; x<xstart+winsize; x++){
      for(int k=0; k<nchannels; k++){
	double diff = im1(x,y,k) - im2(x+dx,y+dy,k);
	sum += diff*diff;
      }
    }
    int x = 0;
    for(x=xstart+hwin; x<xend-hwin-1; x++){
      hpass(x,y) = sum;

      for(int k=0; k<nchannels; k++){
	double diff = im1(x-hwin,y,k) -  im2(x-hwin+dx,y+dy,k);
	sum -= diff*diff;

	diff = im1(x+hwin+1,y,k) -  im2(x+dx+hwin+1,y+dy,k);
	sum += diff*diff;
      }
    }
    hpass(x,y) = sum;
  }

  /*
  Imagef vfilt(1, 2*hwin+1, 1);
  vfilt = 1;
  
  return hpass.convolve(vfilt);
  */

  //Vertical pass.
  for(int x=xstart; x<xend; x++){
    float sum = 0;
    for(int y=ystart; y<ystart+winsize; y++)
      sum += hpass(x,y);
    
    int y = 0;
    for(y=ystart+hwin; y<yend-hwin-1; y++){
      result(x,y) = sum;
      sum += hpass(x,y+hwin+1)-hpass(x,y-hwin);
    }
    result(x,y) = sum;
  }
  return result;  
}

//I0(x,y) = I1(x+dx, y+dy)+Grad(I1)*[dx,dy]
Imagef opticFlowTraditional(const Imagef & im1, const Imagef & im2, const Imagef & estimate, int hwin, Imagef * det){
  Imagef im1gray = (im1.getChannel(0)+im1.getChannel(1)+im1.getChannel(2))*0.3333;
  Imagef im2gray = (im2.getChannel(0)+im2.getChannel(1)+im2.getChannel(2))*0.3333;

  printf("Getting gradient\n");

  Imagef g1 = im1gray.gradient(0);
  Imagef flow(im1.w, im1.h, 2);

  printf("getting copy\n");

  int winsize = 2*hwin+1;

  flow = estimate.copy();

  Matrix A(winsize*winsize, 2);
  Matrix b(winsize*winsize, 1);

  if(det)*det = Imagef(im1.w, im1.h, 1);

  for(int y=hwin; y<im1.h-hwin; y++){
    for(int x=hwin; x<im1.w-hwin; x++){
      int i = 0;
      for(int yy=-hwin; yy<=hwin; yy++){
	for(int xx=-hwin; xx<=hwin; xx++, i++){
	  A(i, 0) = g1(xx+x,yy+y,0);
	  A(i, 1) = g1(xx+x,yy+y,1);	
	}
      }
      Matrix U, V;
      Matrix AtA = (A.transpose()*A);
      (AtA*(1.0/(winsize*winsize))).eigSym(U, V);
      ///double d = AtA.det();
      ///double ratio = std::max(V[0], V[1])/std::min(V[0], V[1]);
      ///printf("det %f %f\n", V[0], V[1]);
      
      if(det)(*det)(x,y) = std::min(V[0], V[1]);

      Matrix Ainv = AtA.inverse()*A.transpose();

      for(int its=0; its<3; its++){
	i = 0;
	double dx = flow(x,y,0);
	double dy = flow(x,y,1);

	for(int yy=-hwin; yy<=hwin; yy++){
	  for(int xx=-hwin; xx<=hwin; xx++, i++){
	    double x2 = std::max(0.0, std::min(dx+x+xx, (double)im2.w-1));
	    double y2 = std::max(0.0, std::min(dy+y+yy, (double)im2.h-1));
	    
	    b[i] = im2gray.bilinear(x2, y2) - im1gray(x+xx, y+yy);
	  }
	}
	Matrix upd = Ainv*b;

	flow(x,y,0) -= upd[0]*0.7;
	flow(x,y,1) -= upd[1]*0.7;
	
	//Matrix res = (A*upd-b);
	//printf("%d,%d  %f\n", x, y, res.dot(res));
      }
    }
  }
  return flow;
}

Imagef opticFlow(const Imagef & im1,const Imagef & im2,int hwin,
		 int x0, int x1, int y0, int y1){
  Imagef bestResult(im1.w, im1.h, 1);
  Imagef bestFrom(im1.w, im1.h, 2);

  bestFrom = 0;
  bestResult = (2*hwin+1)*3+1;

  for(int y=y0; y<y1; y++){
    printf("%d/%d\n", y,y1);
    for(int x=x0; x<x1; x++){
      Imagef result = getSSD_sse(im1, im2, hwin, x, y);
      for(int yy=0; yy<result.h; yy++){
	for(int xx=0; xx<result.w; xx++){
	  if(result(xx,yy)<bestResult(xx,yy)){
	    bestResult(xx,yy) = result(xx,yy);
	    bestFrom(xx,yy,0) = x;
	    bestFrom(xx,yy,1) = y;
	  }
	}
      }
    }
  }
  return bestFrom;
}

Imagef getSSDBrute(const Imagef & im1, const Imagef & im2,
		   int hwin,
		   int dx, int dy){
  int xstart = std::max(0, -dx);
  int xend   = std::min(im1.w, im2.w-dx);

  int ystart = std::max(0, -dy);
  int yend   = std::min(im1.h, im2.h-dy);

  int winsize = 2*hwin+1;
  Imagef result(im1.w, im1.h, 1);
  int nchannels = std::min(3, im1.nchannels);
  
  double maxScore = nchannels*winsize+1;

  result = maxScore;

  for(int y=ystart+hwin; y<yend-hwin; y++){
    for(int x=xstart+hwin; x<xend-hwin; x++){
      double sum = 0;
      for(int yy=-hwin; yy<=hwin; yy++){
	for(int xx=-hwin; xx<=hwin; xx++){
	  for(int k=0; k<nchannels; k++){
	    double diff = im1(x+xx,y+yy,k) - im2(x+xx+dx,y+yy+dy,k);
	    sum+=diff*diff;
	  }
	}
      }
      result(x,y) = sum;
    }
  }
  return result;
}

#ifdef FLOW_MAIN

#include <nmisc/commandline.h>
#include "filters.h"
#include "conversions.h"

int main(int ac, char * av[]){
  double forbackthresh = 10;

  nacb::CommandLine cline("./flow arg1 arg2 arg3");
  cline.registerOption("compose", "Compose two flow fields.");
  cline.registerOption("warp", "Usage --warp disp.flo input.img output.img");
  cline.registerOption("consist", "Consistency operation (check consistency between arg0 arg1 and write results into arg2).  See consist_thresh.");
  cline.registerOption("consist_thresh", "Forward backward consistency threshold", &forbackthresh, 0);
  cline.registerOption("cross_fade", "Cross fade ()");
  cline.registerOption("add_mask", "Add masked flow.");
  cline.parse(ac, av);

  if(cline.gotArgument("compose")){
    if(optind + 2 < ac){
      printf("Composing %s with %s into %s\n", av[optind], av[optind+1], av[optind+2]);
      Imagef U = flow_read_flo(av[optind]);
      Imagef upd = flow_read_flo(av[optind+1]);
      Imagef C = flow_compose(U, upd);
      flow_save_flo(C, av[optind+2]);
    }
    else {
      printf("--compose  Need's two input arguments and an output argument to compose.\n");
      return 1;
    }
    return 0;
  }
  if(cline.gotArgument("add_mask")){
    if(optind + 3 < ac){
      nacb::Imagef d1, d2;
      nacb::Image8 maskImage(av[optind+2]);
      
      d1 = flow_read_flo(av[optind]);
      d2 = flow_read_flo(av[optind+1]);
      
      assert(d1.nchannels == 2 && d2.nchannels == 2);
      assert(d1.w == d2.w && d1.h == d2.h);
      assert(d1.w == maskImage.w && d1.h == maskImage.h);
      
      int use = maskImage.nchannels - 1;
      printf("using mask channel %d\n", use);
      for(int y=0; y<maskImage.h; y++){
	for(int x=0; x<maskImage.w; x++){
	  if(maskImage(x, y, use) > 128){
	    d1(x, y, 0) = d2(x, y, 0);
	    d1(x, y, 1) = d2(x, y, 1);
	  }
	}
      }
      flow_save_flo(d1, av[optind+3]);
    }
    else {
      printf("Need d.flo dover.flo mask.img output.flo");
      return 1;
    }
    return 0;
  }
  if(cline.gotArgument("cross_fade")){
    if(optind + 4 < ac){
      printf("Cross-fading %s %s\n", av[optind], av[optind+1]);
      nacb::Imagef im1(av[optind]), im2(av[optind+1]);
      nacb::Imagef d1, d2;

      d1 = flow_read_flo(av[optind+2]);
      d2 = flow_read_flo(av[optind+3]);

      flow_save_crossfade("/tmp/cf-%04d.tga", 
			  im1, im2, d1, d2, atoi(av[optind+4]));
    }
    else {
      printf("Need the two images, the two flows, and time.\n");
      return 1;
    }
    return 0;
  }

  if(cline.gotArgument("warp")){
    if(optind + 2 < ac){
      printf("Warping %s with %s into %s\n", av[optind], av[optind+1], av[optind+2]);
      Imagef disp = flow_read_flo(av[optind]);
      printf("%dx%dx%d\n", disp.w, disp.h, disp.nchannels);
      Imagef img(av[optind+1]);

      flow_warp_back(disp, img).write(av[optind+2]);
    }
    else {
      printf("--warp  Use disp.flo to warp image.img into output.img");
      return 1;
    }
    return 0;
  }
  if(cline.gotArgument("consist")){
    if(optind + 2 < ac){
      printf("Checking consistency of %s with %s, writing mask into %s\n", av[optind], av[optind+1], av[optind+2]);
      Imagef forward = flow_read_flo(av[optind]);
      Imagef backward = flow_read_flo(av[optind+1]);
      //FIXME: not sure which works better, compositional approach, or basic
      Image8 img = flow_ensure_consistency(forward, backward, forbackthresh, 0);
      
      img.write(av[optind+2]);
    }
    else {
      printf("--warp  Use disp.flo to warp image.img into output.img");
      return 1;
    }
    return 0;
  }

  Imagef im1(av[1]), im2(av[2]);
  
  Imagef forward, backward;

  int hwin = 3;
  
  //im1 = im1.resize(im1.w/2, im1.h/2);
  //im2 = im2.resize(im2.w/2, im2.h/2);
  
  int rgx = 20;
  int rgy = 10;
  
  bool onelevel = true;

  if(onelevel){
    Imagef im1small = im1.resize(im1.w/2, im1.h/2);
    Imagef im2small = im2.resize(im2.w/2, im2.h/2);
    
    forward  = opticFlow(im1small, im2small, hwin, -rgx, rgx, -rgy, rgy);
    backward = opticFlow(im2small, im1small, hwin, -rgx, rgx, -rgy, rgy);
   
    forward = opticFlowTraditional(im1, im2, forward.resize(im1.w, im1.h)*2.0, hwin);
    backward = opticFlowTraditional(im2, im1, backward.resize(im2.w, im2.h)*2.0, hwin);
  }
  else{
    vector<Imagef> im1pyramid;
    vector<Imagef> im2pyramid;
    im1pyramid.push_back(im1);
    im2pyramid.push_back(im2);

    double sc = 0.75;

    for(int its=0; its<3; its++){
      im1pyramid.push_back(im1pyramid[its].resize(im1pyramid[its].w*sc, im1pyramid[its].h*sc));
      im2pyramid.push_back(im2pyramid[its].resize(im2pyramid[its].w*sc, im2pyramid[its].h*sc));    
    }
    
    for(int level = im1pyramid.size()-1; level>=0; level--){
      printf("size %d %d\n", im1pyramid[level].w, im1pyramid[level].h);
      if(level == im1pyramid.size()-1){
	//forward = opticFlowTraditional(im1pyramid[level], im2pyramid[level], im1pyramid[level].copy()*0, hwin);
	//backward = opticFlowTraditional(im2pyramid[level], im1pyramid[level], im2pyramid[level].copy()*0, hwin);
	forward = opticFlow(im1pyramid[level], im2pyramid[level], hwin, -10, 11, -10, 11);
	backward = opticFlow(im2pyramid[level], im1pyramid[level], hwin, -10, 11, -10, 11);
      }
      else {
	forward = forward.resize(im1pyramid[level].w, im1pyramid[level].h)*1.0/sc;
	printf("resizing %d %d\n", forward.w, forward.h);
	forward = opticFlowTraditional(im1pyramid[level], im2pyramid[level], forward, hwin);

	backward = opticFlowTraditional(im2pyramid[level], im1pyramid[level], 
				      backward.resize(im2pyramid[level].w, im2pyramid[level].h)*1.0/sc, hwin);
      }
    }
  }
  ((forward.getChannel(0)-10)*(1.0/20)).write("/tmp/flow-x.ppm");
  
  Imagef warped(im1.w, im1.h, im1.nchannels);
  int badvalue = -400;

  flow_ensure_consistency(forward, backward, forbackthresh, badvalue);

  Imagef xdisp = forward.getChannel(0);
  Imagef ydisp = forward.getChannel(1);

  fill_holes(xdisp, badvalue, 400*1000);
  fill_holes(ydisp, badvalue, 400*1000);
  
  forward.setChannel(0, xdisp);
  forward.setChannel(1, ydisp);


  for(int ti=0; ti<=10; ti++){
    double t = double(ti)/10.0;
    warped = 0;
  
    for(int y=0; y<warped.h; y++){
      for(int x=0; x<warped.w; x++){
	/*int x2 = std::max(0, std::min(x+(int)forward(x,y,0), im2.w-1));
	int y2 = std::max(0, std::min(y+(int)forward(x,y,1), im2.h-1));
	
	int xback = backward(x2,y2,0)+x2;
	int yback = backward(x2,y2,1)+y2;
	
	if(fabs(xback-x)>=2 || fabs(yback-y)>=2){
	  continue;
	}
	*/
	if(forward(x,y,0)==-400 || forward(x,y,1)==-400)continue;
	double xuse = std::max(0.0, std::min(t*forward(x,y,0)+x, (double)im2.w-1));
	double yuse = std::max(0.0, std::min(t*forward(x,y,1)+y, (double)im2.w-1));

	for(int k=0; k<im1.nchannels; k++){
	  warped(x,y,k) = im2.bilinear(xuse, yuse, k);
	}
      }
    }
    
    char fname[1024];
    
    snprintf(fname, 1024, "/tmp/warped-%02d.png", ti);
		     
    warped.write(fname);
  }
  
  Imagef result = forward.getChannel(0)*(warped.getChannel(0)>0);
  double n,x;
  result.getRange(n,x);
  ((result-n)*(1.0/(x-n))).write("/tmp/ssd.png");
  printf("range %f, %f\n", n, x);

  writeRawDepth("/tmp/flow.rfi", forward);

  return 0;
}

#endif
