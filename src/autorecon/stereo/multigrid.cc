#include "multigrid.h"

namespace mg {

  restriction_average::restriction_average(int _w, int _h, 
					   int _neww, int _newh) : 
    w(_w), h(_h), neww(_neww), newh(_newh){
 
    assert(neww < w && newh < h);
  
    locoeff = 0;
    hicoeff = 0;

    construct(w, h, neww, newh);
  }


  restriction_average::~restriction_average(){
    delete [] locoeff; 
    delete [] hicoeff;
  }


  void restriction_average::construct(int _w, int _h, int _neww, int _newh){    
    printf("constructing average %dx%d => %dx%d\n", _w, _h, _neww, _newh);
    
    if(locoeff)delete [] locoeff;
    if(hicoeff)delete [] hicoeff;

    w = _w;
    h = _h;
    neww = _neww;
    newh = _newh;

    double xrat = double(neww)/double(w);
    double yrat = double(newh)/double(h);
    
    double wt = double(neww*newh)/double(w*h);
    
    locoeff = new std::vector<entry>[neww*newh];
    hicoeff = new std::vector<entry>[w*h];

    printf("Creating %d x %d\n", neww, newh);

    for(int yh=0; yh<h; yh++){
      double cy_lo =  double(yh)/double(h);
      double cy_hi = (double(yh)+0.9999)/double(h);
    
      int y_lo = (int)floor(cy_lo*newh), y_hi = (int)floor(cy_hi*newh);
    
      double b = (y_lo==y_hi)?1.0:((y_lo+1)-cy_lo*newh)/yrat;

      for(int xh=0; xh<w; xh++){
	double cx_lo =  double(xh)/double(w);
	double cx_hi = (double(xh)+0.9999)/double(w);
      
	int x_lo = (int)floor(cx_lo*neww), x_hi = (int)floor(cx_hi*neww);
      
	double a = (x_lo==x_hi)?1.0:((x_lo+1)-cx_lo*neww)/xrat;

	//The checks here are so there aren't duplicate entries with zeros.
	if(b>0){
	  if(a>0)
	    locoeff[y_lo*neww + x_lo].push_back(entry(xh, yh, a*b*wt));
	  if(a<1)
	    locoeff[y_lo*neww + x_hi].push_back(entry(xh, yh, (1.0-a)*b*wt));
	}
	if(b<1){
	  if(a>0)
	    locoeff[y_hi*neww + x_lo].push_back(entry(xh, yh, a*(1.0-b)*wt));
	  if(a<1)
	    locoeff[y_hi*neww + x_hi].push_back(entry(xh, yh, (1.0-a)*(1.0-b)*wt));
	}
      }
    }

    //Create transposed version
    for(int i=0; i<neww*newh; i++){
      int xlo = i%neww;
      int ylo = i/neww;

      for(int j=0; j<(int)locoeff[i].size(); j++){
	int y = locoeff[i][j].y;
	int x = locoeff[i][j].x;
	
	hicoeff[y*w + x].push_back(entry(xlo, ylo, locoeff[i][j].val/wt));
      }
    }

    print();
  }


  nacb::Imagef restriction_average::prolong(const nacb::Imagef & input){
    assert(input.w == neww && input.h == newh);
   
    nacb::Imagef im(w, h, input.nchannels);

    for(int y=0; y<h; y++){
      for(int x=0; x<w; x++){
	float sum[input.nchannels];
	memset(sum, 0, sizeof(float)*input.nchannels);

	int ind = y*w + x;
	for(int i=0; i<(int)hicoeff[ind].size(); i++){

	  for(int k=0; k<input.nchannels; k++){
	    sum[k] += input(hicoeff[ind][i].x, hicoeff[ind][i].y, k) * hicoeff[x][i].val;

	  }
	}
	memcpy(&im(x, y, 0), sum, sizeof(float)*input.nchannels);
      }
    }
    return im;
  }

  nacb::Imagef restriction_average::apply(const nacb::Imagef & input){
    nacb::Imagef im(neww, newh, input.nchannels);

    for(int y=0; y<newh; y++){
      for(int x=0; x<neww; x++){
	float sum[input.nchannels];
	memset(sum, 0, sizeof(float)*input.nchannels);

	int ind = y*neww + x;
	for(int i=0; i<(int)locoeff[ind].size(); i++){
	  for(int k=0; k<input.nchannels; k++){
	    sum[k] += input(locoeff[ind][i].x, locoeff[ind][i].y, k) * locoeff[x][i].val;
	  }
	}
	memcpy(&im(x, y, 0), sum, sizeof(float)*input.nchannels);
      }
    }
    return im;
  }


  void restriction_average::print(){
    for(int x=0; x<neww; x++){
      printf("%d from %d\n", x, (int)locoeff[x].size());
      for(int i=0; i<(int)locoeff[x].size(); i++){
	printf("%dx%d (%lf) ", locoeff[x][i].x, locoeff[x][i].y, locoeff[x][i].val);
      }
      printf("\n");
    }
  }


};
