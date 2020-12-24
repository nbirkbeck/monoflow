#include <stdio.h>
#include <stdio.h>
#include <nimage/image.h>
#include <nmath/vec2.h>
#include <nmath/vec3.h>
#include <nmath/matrix.h>
#include <nmath/sparsematrix.h>
#include <vector>
#include <queue>
#include <assert.h>

using namespace std;
using namespace nacb;

Imagef fast_median_reject_pass1(const Imagef & im,
			  int hs, 
			  int x1, int x2,
			  float reject){
  int ndisp = (x2-x1)+1;
  int * bins = new int[ndisp];
  memset(bins, 0, sizeof(int)*ndisp);
  
  int pivot = 0;
  int winsize = 2*hs+1;
  int cnt = 0;
  Imagef result(im.w, im.h, 1);

  for(int y=hs; y<im.h-hs; y++){
    memset(bins, 0, sizeof(int)*ndisp);
    cnt = 0;

    for(int yy=-hs; yy<=hs; yy++){
      for(int xx=0; xx<winsize; xx++){
	if(im(xx,y+yy)!=reject){
	  bins[int(im(xx,y+yy)-x1)]++;
	  cnt++;
	}
      }
    }
    int sum = 0;
    pivot = 0;
    
    while(sum<cnt/2 && pivot<ndisp){
      sum += bins[pivot++];
    }
    //x1's value is the median
    //median(hs,y) <- pivot+x1
    for(int x=0; x<hs; x++)
      result(x, y) = pivot+x1-1;

    result(hs, y) = pivot+x1-1;

    for(int x=hs+1; x<im.w-hs; x++){
      int offs = 0;
      for(int yy=-hs; yy<=hs; yy++){

	if(im(x-hs-1, y+yy)!=reject){
	  if(im(x-hs-1, y+yy)<pivot+x1)
	    offs++;
	  bins[int(im(x-hs-1, y+yy)-x1)]--;
	  cnt--;
	}
	if(im(x+hs, y+yy)!=reject){
	  if(im(x+hs, y+yy)<pivot+x1)
	    offs--;

	  bins[int(im(x+hs, y+yy)-x1)]++;
	  cnt++;
	}
      }
      sum -= offs;
      
      //if(offs<0){
	while(pivot>0  && sum>=cnt/2){
	  int newsum = sum - bins[pivot-1];
	  sum = newsum;
	  pivot--;
	}
	if(cnt==0)
	  pivot=0;
	//}
      while(sum<cnt/2 && pivot<ndisp)
	sum += bins[pivot++];

      
      /*     
      int sum2 = 0;
      int total = 0;
      int pivot2 = 0;
      for(int k=0; k<pivot; k++){
	sum2 += bins[k];
      }
      for(int k=0; k<ndisp; k++){
	total += bins[k];
	if(total<cnt/2)
	  pivot2++;
      }
      if(x<100 && (sum2!=sum || total!=cnt || pivot !=pivot2))
	printf("sum %d, %d   (%d,%d)   %d,%d [%d,%d,%d]\n", sum, sum2, cnt, total, pivot, pivot2, bins[pivot-1], bins[pivot], bins[pivot+1]);
      */
      result(x, y) = pivot+x1-1;
      
    }

    for(int x=im.w-hs; x<im.w; x++)
      result(x, y) = pivot+x1-1;
  }
  
  for(int y=0; y<hs; y++){
    for(int x=0; x<im.w; x++){
      result(x,y) = im(x,y);
      result(x,im.h-y-1) = im(x,im.h-y-1);
    }
  }

  delete [] bins;
  return result;
}


Imagef fast_median_reject(const Imagef & im,
			  int hs, 
			  int x1, int x2,
			  float reject){
  Imagef med = fast_median_reject_pass1(im, hs, x1, x2, reject);
  vector<float> diffs(im.w*im.h, 10e10);
  int cnt = 0;
  for(int i=0; i<im.w*im.h; i++)
    if(im.data[i]!=reject)
      diffs[cnt++] = (fabs(im.data[i]-med.data[i]));
  
  std::sort(diffs.begin(), diffs.end());
  double thresh = diffs[cnt/2];

  int rej = 0;
  for(int y=0; y<im.h; y++){
    for(int x=0; x<im.w; x++){
      if(fabs(im(x, y)-med(x,y))>2.0*thresh){
	med(x,y) = reject;
	rej++;
      }
      else {
	med(x,y) = im(x,y);
      }
    }
  }
  printf("rejected %d %f %f\n", rej, thresh, thresh);
  return med;
}


Imagef median_reject(const Imagef & im, 
		     int hs, float reject){
  Imagef result = im.copy();

  for(int y=hs; y<im.h-hs; y++){
    for(int x=hs; x<im.w-hs; x++){
      vector<float> vals;
      for(int yy=-hs; yy<=hs; yy++){
	for(int xx=-hs; xx<=hs; xx++){
	  if(im(x+xx,y+yy)!=reject)
	    vals.push_back(im(x+xx,y+yy));
	}
      }
      std::sort(vals.begin(), vals.end());

      if(vals.size()>2){
	float med = vals[vals.size()/2];
	
	for(int i=0; i<(int)vals.size(); i++){
	  vals[i] = fabs(vals[i] - med);
	}
	std::sort(vals.begin(), vals.end());

	float thresh = std::max(1.6*vals[vals.size()/2],1.0);
	if(fabs(im(x,y)-med)>thresh){
	  result(x,y) = reject;
	}
      }
    }
  }
  return result;
}

Imagef trilateral(const Imagef & im, 
		  const Imagef & sc,
		  int hs, float reject){
  Imagef result = im.copy();

  Imagef gauss(2*hs+1, 2*hs+1, 1);
  double sum = 1;
  
  {
    gauss(hs, hs) = 1.0;    
    for(int y=0; y<=hs; y++){
      for(int x=0; x<=hs; x++){
	if (x==0 && y==0) continue;
	double val = exp(-x*x/32.0-y*y/32);
	gauss( x+hs, y+hs) = val;
	gauss(-x+hs, y+hs) = val;
	gauss(-x+hs,-y+hs) = val;
	gauss( x+hs,-y+hs) = val;
	sum +=4*val;
      }
    }
    gauss *= (1.0/sum);
  }

  for(int y=hs; y<im.h-hs-1; y++){
    for(int x=hs; x<im.w-hs-1; x++){
      double w = 0;
      double c = im(x,y);
      double sum = 0;

      if(im(x,y)==reject)continue;

      for(int xx=-hs; xx<=hs; xx++){
	for(int yy=-hs; yy<=hs; yy++){
	  if(im(x+xx,y+yy) != reject){
	    double d = im(x+xx,y+yy)-c;
	    double r = exp(-d*d/(2.0*10.0*10.0)); // was /4/4
	    double wt = r*gauss(hs+xx, hs+yy)*std::max(std::min(sc(x+xx,y+yy), 1.0f), 0.0f);
	    sum += wt*float(im(x+xx,y+yy));
	    w += wt;
	  }
	}
      }
      if(w>0)
	sum*=(1.0/w);
      else sum = reject;
      result(x,y,0)=(sum);
    }
  }
  return result;
}

void fill_hole(Imagef & disp, Image8 & mask, int bx, int by, int bw, int bh, int hole){
  //Solve laplacian s.t boundary conditions
  int n = 0;
  Imagef var(bw,bh,1);

  for(int yy=by; yy<by+bh; yy++){
    for(int xx=bx; xx<bx+bw; xx++){
      var(xx-bx, yy-by) = n;
      n += mask(xx,yy)==hole;
    }
  }

  Matrix b = Matrix::zeros(n,1);
  vector<vector<pair<int,double> > > cols(n);


  for(int yy=by; yy<by+bh; yy++){
    for(int xx=bx; xx<bx+bw; xx++){
      if(mask(xx,yy)!=hole)continue;

      int xm = xx-1, xp = xx+1, ym = yy-1, yp = yy+1;
      int i = (int)var(xx-bx, yy-by);
      double sum = 0;

      if(xm<0)
	;//ignoring for now
      else if(mask(xm,yy)==hole){
	int j = (int)var(xm-bx, yy-by);
	cols[j].push_back(pair<int,double>(i, 1.0));
	sum++;
      }
      else {
	b[i] -= disp(xm,yy);
	sum++;
      }

      if(ym<0)
	;//ignoring for now
      else if(mask(xx,ym)==hole){
	int j = (int)var(xx-bx, ym-by);
	cols[j].push_back(pair<int,double>(i, 1.0));
	sum++;
      }
      else {
	b[i] -= disp(xx,ym);
	sum++;
      }

      if(xp>=disp.w)
	;//ignoring for now
      else if(mask(xp,yy)==hole){
	int j = (int)var(xp-bx, yy-by);
	cols[j].push_back(pair<int,double>(i, 1.0));
	sum++;
      }
      else {
	b[i] -= disp(xp,yy);
	sum++;
      }

      if(yp>=disp.h)
	;//ignoring for now
      else if(mask(xx,yp)==hole){
	int j = (int)var(xx-bx, yp-by);
	cols[j].push_back(pair<int,double>(i, 1.0));
	sum++;
      }
      else {
	b[i] -= disp(xx,yp);
	sum++;
      }
      
      cols[i].push_back(pair<int,double>(i, -sum));
    }
  }
  

  SparseMatrix A(n, cols);
  Matrix x = SparseMatrix::umfsolve(A, b);

  //x.printMatlab("x");

  int ind = 0;
  for(int yy=by; yy<by+bh; yy++)
    for(int xx=bx; xx<bx+bw; xx++)
      if(mask(xx,yy) == hole)
	disp(xx,yy) = x[ind++];
}

void fill_holes(nacb::Imagef & disp, 
		float reject, int maxsize){
  const int inq = 1;
  const int current_hole = 3;  
  const int processed_hole = 2;

  Image8 visited(disp.w, disp.h, 1);
  visited = 0;
  
  for(int y=1; y<disp.h-1; y++){
    for(int x=1; x<disp.w-1; x++){
      if(disp(x,y) == reject && !visited(x,y)){
	std::queue<Vec2<int> > q;
	vector<Vec2<int> > points;

	q.push(Vec2<int>(x,y));
	
	//Could also find boundary here instead
	while(!q.empty()){
	  Vec2<int> pos = q.front();
	  q.pop();
	  int px = int(pos.x);
	  int py = int(pos.y);

	  if(disp(px, py) != reject){
	    printf("bad stuff in the queue\n");
	    printf("q.size()=%ld\n", q.size());
	  }

	  visited(px,py) = current_hole; //!< Currently processing portion
	  points.push_back(pos);
	  
	  if((int)points.size()>maxsize)break;

	  if(px+1<disp.w && !visited(px+1, py) && disp(px+1, py)==reject){
	    q.push(Vec2<int>(px+1, py));
	    visited(px+1, py) = inq;
	  }
	  
	  if(py+1<disp.h && !visited(px, py+1) && disp(px, py+1)==reject){
	    q.push(Vec2<int>(px, py+1));
	    visited(px, py+1) = inq;
	  }

	  if(px>0 && !visited(px-1, py) && disp(px-1, py)==reject){
	    q.push(Vec2<int>(px-1, py));
	    visited(px-1, py) = inq;
	  }

	  if(py>0 && !visited(px, py-1) && disp(px, py-1)==reject){
	    q.push(Vec2<int>(px, py-1));
	    visited(px, py-1) = inq;
	  }
	}

	//In this case, we found a hole...try and see if it can be filled
	if((int)points.size()<=maxsize){
	  Vec2<int> mn(disp.w, disp.h), mx(0,0);
	  for(int i=0; i<(int)points.size(); i++){
	    mn = mn.min(points[i]);
	    mx = mx.max(points[i]);
	  }
	  int bx = (int)mn.x;
	  int by = (int)mn.y;
	  
	  int bw = (int)mx.x - bx + 1;
	  int bh = (int)mx.y - by + 1;

	  //printf("filling hole at %d,%d to %d,%d   size=%d\n", bx, by, bw, bh, points.size());

	  vector<Vec3f> boundary_vals;
	  //try to find boundary
	  for(int yy=-1; yy<=bh; yy++){
	    int y2 = yy + by;
	    if(y2<0 || y2>=disp.h)continue;

	    for(int xx=-1; xx<=bw; xx++){
	      int x2 = xx + bx;
	      if(x2<0 || x2>=disp.w)continue;
	      
	      //Find some non-hole pixel
	      if(!visited(x2, y2) && disp(x2, y2)!=reject){
		bool border = false;
		
		for(int i=0; i<4; i++){
		  int ox = x2+((i<=1 && i&1)?1:-1);
		  int oy = y2+((i>=2 && i&0x2)?1:-1);
		  
		  if(ox>=0 && ox<disp.w && oy>=0 && oy<disp.h 
		     && visited(ox, oy) == current_hole)
		    border = true;
		}
		
		if(border)
		  boundary_vals.push_back(Vec3f(x2, y2, disp(x2, y2)));
	      }
	    }
	  }
	  if(boundary_vals.size()>=2){

	    Vec3f mean(0,0,0);
	    for(int i=0; i<(int)boundary_vals.size(); i++) mean += boundary_vals[i];
	    mean *= 1.0/boundary_vals.size();

	    Matrix covar(3, boundary_vals.size());
	    for(int i=0; i<(int)boundary_vals.size(); i++){
	      Vec3f diff = boundary_vals[i] - mean;
	      covar(0, i) = diff.x;
	      covar(1, i) = diff.y;
	      covar(2, i) = diff.z;
	    }

	    ///covar.printMatlab("covar");

	    covar = covar*covar.transpose();
	    covar *= 1.0/boundary_vals.size();

	    ///covar.printMatlab("covar");

	    Matrix U,S,V;
	    covar.Lsvd(U,S,V);

	    Vec3f n(U(0,2), U(1,2), U(2,2));
	    double avgdist = 0;
	    for(int i=0; i<(int)boundary_vals.size(); i++){
	      double dist = n.dot(boundary_vals[i]-mean);
	      avgdist += dist*dist;
	    }
	    avgdist = sqrt(avgdist/boundary_vals.size());

	    ///printf("avgdist %f\n", avgdist);

	    //Can we fill it? yes we can. FIXME: threshold on avgdist
	    if(avgdist < 2.0 || n.len()<1e-5) {
	      fill_hole(disp, visited, bx, by, bw, bh, current_hole);
	    }
	  }
	}
	//Set all points in the hole to be processed
	for(int i=0; i<(int)points.size(); i++){
	  visited(points[i].x, points[i].y) = processed_hole;
	}
      }
    }
  }
}
