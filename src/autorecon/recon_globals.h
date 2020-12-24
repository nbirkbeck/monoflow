#ifndef RECON_GLOBALS_H
#define RECON_GLOBALS_H
#include <nmath/matrix.h>
#include <nimage/image.h>
#include <nmath/vec3.h>
#include <nmath/vec2.h>
#include <vector>



nacb::Vec2d distort(const nacb::Matrix & k,double a,double b);
nacb::Vec2d intrinsic(const nacb::Matrix & K,
		      const nacb::Vec2d & p);

template<class T>
nacb::Image<T> rectify(nacb::Image<T> & img, nacb::Matrix & K,nacb::Matrix & kc);


template <class T>
nacb::Vec2<T> project(const nacb::Matrix & P, const nacb::Vec3<T> & pt, T * zout = 0);

nacb::Matrix projectPoints(const nacb::Matrix & P,const nacb::Matrix & X);

template <class T>
nacb::Vec3<T> backProject(const nacb::Matrix & KRinv, const nacb::Matrix & Kt, T x, T y, T z);

template <class T>
nacb::Vec3<T> backProject(const nacb::Matrix & P, T x, T y, T z){
  return backProject(P.submatrix(0, 0, 3, 3).inverse(),
		     P.submatrix(0, 3, 3, 1), x, y, z);
}


double cubeeval(double a,double b,double c,double d,double x);
int cuberoots(double a,double b,double c,double d,
	      double & x1,double & x2,double & x3);
void unique_set(int * cset,int max,int n);
std::vector<int> unique(std::vector<int> & input);

nacb::Image8 rgb2gray(nacb::Image8 & im);

void gradient(const nacb::Image8 & gim,nacb::Imagef & dx,nacb::Imagef & dy);
void gauss_sep(double sigma,nacb::Imagef & gx,nacb::Imagef & gy);

std::vector<nacb::Vec2f > harris(nacb::Image8 & im,int maxpts,int subpix=0);

//From svn/thesis/calibration, in recon_planar_calibrate
nacb::Matrix getExtr(const nacb::Matrix & R,const nacb::Matrix & t);
nacb::Matrix projectPoints(const nacb::Matrix & A,
			   const nacb::Matrix & dist,
			   const nacb::Matrix & extr,
			   const nacb::Matrix & pts);

nacb::Matrix normalizeProjectedPoints(const nacb::Matrix & Ainv,const nacb::Matrix & pts);

nacb::Matrix normalizeProjectedPoints(const nacb::Matrix & Ainv, 
				      const nacb::Matrix & d,
				      const nacb::Matrix & pts);

void factorProjectionMatrix(const nacb::Matrix & P, nacb::Matrix & K, nacb::Matrix & E);


class MatchAndScore{
public:
  int index1;
  int index2;
  double score;
  MatchAndScore(int i1=0,int i2=0,double sc=0){
    index1=i1;
    index2=i2;
    score=sc;
  }  
  int operator<(const MatchAndScore & pad) const{
    return score<pad.score;
  }
  int operator<=(const MatchAndScore & pad) const{
    return score<=pad.score;
  }
  int operator>(const MatchAndScore & pad) const{
    return score>pad.score;
  }
  int operator>=(const MatchAndScore & pad) const{
    return score>=pad.score;
  }
  int operator==(const MatchAndScore & pad) const{
    return score==pad.score;
  }
};


struct IndexAndScore{
  int index;
  double score;
  
  int operator<(const IndexAndScore & pad) const{
    return score<pad.score;
  }
  int operator<=(const IndexAndScore & pad) const{
    return score<=pad.score;
  }
  int operator>(const IndexAndScore & pad) const{
    return score>pad.score;
  }
  int operator>=(const IndexAndScore & pad) const{
    return score>=pad.score;
  }
  int operator==(const IndexAndScore & pad) const{
    return score==pad.score;
  }
};

struct PointAndScore{
  nacb::Vec2f point;
  double score;
  
  int operator<(const PointAndScore & pad) const{
    return score<pad.score;
  }
  int operator<=(const PointAndScore & pad) const{
    return score<=pad.score;
  }
  int operator>(const PointAndScore & pad) const{
    return score>pad.score;
  }
  int operator>=(const PointAndScore & pad) const{
    return score>=pad.score;
  }
  int operator==(const PointAndScore & pad) const{
    return score==pad.score;
  }
};

#endif
