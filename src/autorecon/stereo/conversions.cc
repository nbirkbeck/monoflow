#include "conversions.h"
#include <nmath/vec3.h>
#include <nmath/vec4.h>
#include "../recon_globals.h"
using namespace nacb;


void disparityToDepthHW(const Matrix & Pn1,
			const Matrix & Pn2,
			double minx1, double minx2,
			double disp_reject, 
			double depth_reject){
  fprintf(stderr, "Not written yet..\n");
}

//f*(x/z)+d = f*(x-baseline)/(z)
// ==> d = -f*baseline/z ==> z = -f*baseline/d
// ==> d = -f*baseline/z + (minx2-minx1)
// ==> f*baseline =  z*(minx2-minx1)-d*z
Imagef disparityToDepth(const Imagef & disp,
			const Matrix & Pn1,
			const Matrix & Pn2,
			double minx1, double minx2,
			double disp_reject,
			double depth_reject){
  
  Matrix K1, E1, K2, E2;
  factorProjectionMatrix(Pn1, K1, E1);
  factorProjectionMatrix(Pn2, K2, E2);

  Matrix cc2 = (E2.inverse()).getColumn(3);
  Matrix diff = E1*cc2;
  double baseline = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
  if(diff[0]<0)
    baseline *= -1;

  double f = (K1(0,0) + K2(0,0))/2.0;
  Imagef depth(disp.w, disp.h, 1);

  for(int y=0; y<disp.h; y++){
    for(int x=0; x<disp.w; x++){
      double dp = disp(x,y);
      if(dp == disp_reject)
	depth(x,y) = depth_reject;
      else {
	double d = dp + minx2 - minx1;
	depth(x,y) = -f*baseline/d;
      }
    }
  }
  //printf("baseline: %f\n", baseline);
  //E1.printMatlab("E1");
  //E2.printMatlab("E2");
  //diff.printMatlab("diff");
  return depth;
}

Imagef warpDepth(int w, int h,
		 Imagef & disp,
		 Matrix & E,
		 Matrix & Pnew,
		 Matrix & H,
		 double minx, double miny){
  Matrix sinv = Pnew.submatrix(0,0,3,3).inverse();
  Matrix sub = Pnew.submatrix(0,3,3,1);
  Matrix A = Matrix::eye(4,4);

  A.setSubmatrix(0,0, sinv);

  Matrix depthTransform = E.getRow(2)*A*Matrix::trans(-sub[0], -sub[1], -sub[2]);
  Vec4d dotter(depthTransform[0], depthTransform[1], depthTransform[2], depthTransform[3]);

  Imagef depth(w,h,1);

  for(int y=0; y<h; y++){
    for(int x=0; x<w; x++){
      double hx = H(0,0)*x + H(0,1)*y + H(0,2);
      double hy = H(1,0)*x + H(1,1)*y + H(1,2);
      double hz = H(2,0)*x + H(2,1)*y + H(2,2);

      hx = hx/hz;
      hy = hy/hz;

      double z = disp.bilinear(hx-minx, hy-miny, 0);

      /*
      Matrix t(3,1);
      
      t[0] = hx*z - sub[0];
      t[1] = hy*z - sub[1];
      t[2] = z - sub[2];
      
      t = sinv*t;
      Matrix t4(4,1);
      t4.setSubmatrix(0,0,t);
      t4[3] = 1;
      
      t4 = Pnew * t4;
      t4[0] /= t4[2];
      t4[1] /= t4[2];
      printf("%f %f %f   %f %f %f\n", t4[0], t4[1], t4[2], hx, hy, z);
      */
      depth(x,y) = dotter.x*hx*z + dotter.y*hy*z + dotter.z*z + dotter.w;
    }
  }
  return depth;
}

Imagef depthToPoints(const Imagef & depth,
		     const Matrix & P,
		     double minx, double miny){
  Matrix A = P.submatrix(0,0,3,3).inverse();
  Matrix b = P.submatrix(0,3,3,1);

  Vec3d r1(A(0,0), A(0,1), A(0,2));
  Vec3d r2(A(1,0), A(1,1), A(1,2));
  Vec3d r3(A(2,0), A(2,1), A(2,2));

  Imagef points(depth.w, depth.h, 3);

  for(int y=0; y<depth.h; y++){
    for(int x=0; x<depth.w; x++){
      double z = depth(x, y);
      Vec3d p = Vec3d((x+minx)*z-b[0], (y+miny)*z-b[1], z-b[2]);

      points(x,y,0) = r1.dot(p);
      points(x,y,1) = r2.dot(p);
      points(x,y,2) = r3.dot(p);
    }
  }
  return points;
}

void savePoints(const Imagef & depth,
		const Imagef & points,
		const char * fname,
		double depthReject,
		double thresh){
  FILE * file = fopen(fname, "w");

  for(int y=0; y<points.h; y++){
    for(int x=0; x<points.w; x++){
      fprintf(file, "v %f %f %f\n", points(x, y, 0),
	      points(x, y, 1), points(x, y, 2));
    }  
  }
  
  for(int y=0; y<points.h-1; y++){
    for(int x=0; x<points.w-1; x++){
      if(depth(x,y) == depthReject ||
	 depth(x+1,y) == depthReject ||
	 depth(x,y+1) == depthReject)continue;
            
      if(fabs(depth(x,y) - depth(x+1,y))<thresh &&
	 fabs(depth(x,y) - depth(x+1,y+1))<thresh &&
	 fabs(depth(x,y) - depth(x,y+1))<thresh){
	fprintf(file, "f %d %d %d\n", 
		y*points.w + x  +1,
		(y+1)*points.w + x  +1,
		(y+1)*points.w + x+1  +1);
	fprintf(file, "f %d %d %d\n", 
		y*points.w + x  +1,
		(y+1)*points.w + x+1  +1,
		y*points.w + x+1  +1);
      }
    }  
  }
  fclose(file);
}

void writeRawDepth(const char * imname, const Imagef & im){
  FILE * file = fopen(imname, "wb");
  uint32_t dims[4] = {(uint32_t)im.w, (uint32_t)im.h, (uint32_t)im.nchannels,0};

  fwrite(dims, sizeof(uint32_t), 4, file);
  fwrite(im.data, sizeof(float), im.w*im.h*im.nchannels, file);
  
  fclose(file);
}

Imagef readRawDepth(const char * imname){
  uint32_t dims[4];

  FILE * file = fopen(imname, "rb");
  Imagef im(0,0,0);

  if(!file)
    return im;

  if(4 != fread(dims, sizeof(uint32_t), 4, file))
    return im;
  
  im = Imagef(dims[0], dims[1], dims[2]);
  const size_t r = fread(im.data, sizeof(float), im.w*im.h*im.nchannels, file);
  if ((int)r != im.w*im.h*im.nchannels) {
    fprintf(stderr, "%s:%d: Read %ld, expected %d\n",
            __FILE__, __LINE__, r, im.w*im.h*im.nchannels);
  }
  
  fclose(file);
  return im;
}

#ifdef CONVERSION_MAIN_H

#include <getopt.h>
#include "utils.h"
#include "../clbfile.h"
#include "zncc_stereo.h"
#include "filters.h"

int main(int ac, char * av[]){
  Image8 images[2];
  Matrix K1(3,3), E1(4,4), K2(3,3), E2(4,4);
  Matrix T1(3,3), T2(3,3), Pn1(3,4), Pn2(3,4);
  Matrix d1(5,1), d2(5,1);
  double minx1, minx2, minx, maxx, miny, maxy;
  double maxx1, maxx2;
  int hwin = 2;
  bool downsample = false;
  double minz = 0.5;
  double maxz = 3.5;
  bool verbose = false;
  bool subpixel = false;

  int longind, ind;
  struct option opts[] = {{"minz",1,0,'n'},
			  {"maxz",1,0,'x'},
			  {"downsample",0,0,'d'},
			  {"hwin",1,0,'h'},
			  {"verbose",0,0,'v'},
			  {"subpixel",0,0,'s'},
			  {0,0,0,0}};
			  

  while((ind = getopt_long(ac, av, "n:x:dh:vs", opts, &longind)) >= 0){
    switch(ind){
    case 's':
      subpixel = true;
      break;
    case 'h':
      hwin = atoi(optarg);
      break;
    case 'n':
      minz = atof(optarg);
      break;
    case 'x':
      maxz = atof(optarg);
      break;
    case 'd':
      downsample = true;
      break;
    case 'v':
      verbose = true;
      break;
    default:
      printf("unknown option %c (int:%d)\n", ind, ind);
      break;
    }
  }

  if(ac - optind  != 2){
    printf("Need the two images as input arguments\n\n");
    printf("This program generates a depth mesh from the two input images\n");
    printf("using a fixed depth range\n");
    return 1;
  }

  if(verbose){
    printf("Starting conversion sample program.  Settings:\n");
    printf("\tinput images: %s, %s\n", av[optind], av[optind+1]);
    printf("\thwin: %d, minz: %f, maxz: %f\n", hwin, minz, maxz);
    printf("\tdownsample: %d\n", downsample);
    printf("\n");
  }
 
  images[0] = Image8(av[optind]);
  images[1] = Image8(av[optind+1]);

    
  if(!ClbFile::read((std::string(av[optind])+".clb").c_str(),
                    K1.data, E1.data, Pn1.data, d1.data)){
    fprintf(stderr, "Error reading calibration file for %s\n", av[optind]);
    exit(2);
  }
  if(!ClbFile::read((std::string(av[optind+1])+".clb").c_str(), 
                    K2.data, E2.data, Pn2.data, d2.data)){
    fprintf(stderr, "Error reading calibration file for %s\n", av[optind+1]);
    exit(2);
  }
  
  double sc = downsample ? 0.5 : 1;
  Matrix sc_mat = Matrix::eye(3,3);
  sc_mat(0,0) = sc_mat(1,1) = sc;
  
  K1 = sc_mat*K1;
  K2 = sc_mat*K2;
  
  images[0] = images[0].resize(sc*images[0].w, sc*images[0].h);
  images[1] = images[1].resize(sc*images[1].w, sc*images[1].h);
  
  rectify(K1, E1, K2, E2, T1, T2, Pn1, Pn2);
  findBounds(images[0], T1, images[1], T2, minx, maxx, miny, maxy);
  minx1 = minx2 = 10000;
  maxx1 = maxx2 = 0;
  findBounds(images[0], T1, minx1, maxx1);
  findBounds(images[1], T2, minx2, maxx2);

  if(verbose){
    printf("Bounds: %f, %f\n", minx1, maxx1);
    printf("Bounds: %f, %f\n", minx2, maxx2);
  }
  Image8 rect0 = applyHomography(images[0], T1, minx1, maxx1, miny, maxy);
  Image8 rect1 = applyHomography(images[1], T2, minx2, maxx2, miny, maxy);

  Imagef r0f, r1f;
  r0f = rect0;
  r1f = rect1;

  Vec2d drange = getDisparityRange(Pn1, Pn2, minx1, minx2, minz, maxz);
  int x1 = (int)floor(drange.x);
  int x2 = (int)ceil(drange.y);

  if(verbose)
    printf("x1:%d x2:%d\n", x1, x2);
  
  Imagef score_left;
  Imagef disp_left = getMultiScaleDisparity(r0f, r1f, hwin, x1, x2, &score_left);

  Imagef score_right;
  Imagef disp_right = getMultiScaleDisparity(r1f, r0f, hwin, -x2, -x1, &score_right);

  disp_left = median_reject(disp_left, 7, x1-1.0);
  disp_right = median_reject(disp_right, 7, -x2-1.0);

  ensureForwardBackward(disp_left, disp_right, 3.0, x1-1);
  ensureForwardBackward(disp_right, disp_left, 3.0, -x2-1);
  
  if(subpixel){
    if(verbose)
      printf("Doing subpixel.\n");

    zncc_subpixel(r0f, r1f, disp_left, hwin, x1-1, 3, true);
    zncc_subpixel(r1f, r0f, disp_right, hwin, -x2-1, 3, true);
  }
  
  disp_left = median_reject(disp_left, 7, x1-1.0);
  fill_holes(disp_left, x1-1, 100*100);
  disp_left = trilateral(disp_left, score_left, 7, x1-1);

  disp_right = median_reject(disp_right, 7, -x2-1);
  fill_holes(disp_right, -x2-1, 100*100);
  disp_right = trilateral(disp_right, score_right, 7, -x2-1);

  Imagef depth_left = disparityToDepth(disp_left, Pn1, Pn2, minx1, minx2, x1-1, 0.0);
  Imagef depth_right = disparityToDepth(disp_right, Pn2, Pn1, minx2, minx1, -x2-1, 0.0);
  /*
  Matrix T1tr = Matrix::eye(3,3);
  
  T1tr(0,2) = -minx1;
  T1tr(1,2) = -miny;
  
  Matrix T1_inv = (T1tr*T1).inverse();
  Imagef depth = applyHomography(depth_warped, T1_inv, 0.0, (double)images[0].w, 0, (double)images[0].h);
  */
  depth_left = warpDepth(images[0].w, images[0].h, depth_left, E1, Pn1, T1, minx1, miny);
  depth_right = warpDepth(images[1].w, images[1].h, depth_right, E2, Pn2, T2, minx2, miny);

  if(verbose){
    printf("writing output images\n");
  }

  ((depth_left-minz)*1.0/(maxz-minz)).write("/tmp/depth_left.png");
  ((depth_right-minz)*1.0/(maxz-minz)).write("/tmp/depth_right.png");


  writeRawDepth("/tmp/depth_left.rfi", depth_left);
  writeRawDepth("/tmp/depth_right.rfi", depth_right);


  Matrix P1 = K1*Matrix::eye(3,4)*E1;
  Imagef points_left = depthToPoints(depth_left, P1, 0.0, 0.0);// minx1, miny);

  Matrix P2 = K2*Matrix::eye(3,4)*E2;
  Imagef points_right = depthToPoints(depth_right, P2, 0.0, 0.0);//minx2, miny);

  if(verbose){
    printf("writing objects\n");
  }

  savePoints(depth_left, points_left, "/tmp/left.obj", 0.0, (maxz-minz)/60.0);
  savePoints(depth_right, points_right, "/tmp/right.obj", 0.0, (maxz-minz)/60.0);

  return 0;
}


#endif
