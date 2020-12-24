#ifndef IK_TRIANGLE_H
#define IK_TRIANGLE_H

#include <nmath/vec3.h>
#include <nmath/matrix.h>

class Triangle{
public:
  nacb::Vec3f v1, v2, v3, n;
  nacb::Vec3f i1, i2;     //!< Cached direction to project when checking for ray/triangle intersection.
  nacb::Vec3f e1, e2, e3; //!< Edges (v2-v1), (v3-v1), and, strangely (v3-v2) 
  double d; //! Cached v1.dot(n), e.g., n.dot(x) = d  for points on plane

  ///Invalid triangle
  Triangle(){}

  /*
    Construct the triangle and cache stuff for bary-centric projection onto the triangle.
  */
  Triangle(const nacb::Vec3f & _v1,
	   const nacb::Vec3f & _v2,
	   const nacb::Vec3f & _v3) : v1(_v1), v2(_v2), v3(_v3){
    e1 =(v2-v1);
    e2 =(v3-v1);
    e3 =(v3-v2);

    n = e1.cross(e2);
    n.normalize();
    
    d = v1.dot(n);

    double A = e1.dot(e1);
    double B = e1.dot(e2);
    double C = e2.dot(e2);

    double det = A*C - B*B;
    
    i1 = (e1*C - e2*B)*(1.0/det);
    i2 = (e2*A - e1*B)*(1.0/det);
  }

  double set_less(double  dist, double test, const nacb::Vec3d & wts, nacb::Vec3d * weights) const {
    if(test<dist){
      if(weights)*weights = wts;
      return test;
    }
    return dist;
  }

  //FIXME: This is a ray function, not a triangle function.
  static double distanceToRay(const nacb::Vec3f & o, 
			      const nacb::Vec3f & k,
			      const nacb::Vec3f & pt){
    return (pt - (k*((pt-o).dot(k)) + o)).len();
  }

  /** 
      Compute the distance to the suplied ray.  Doesn't check for intersections.
      Returns the minimum distance and sets the vweights for the vertices (bary-centric)
      if vweights is not null.
  */
  double distanceToRay(const nacb::Vec3f & o,
		       const nacb::Vec3f & k,
		       nacb::Vec3d * vweights = 0){
    double dist = 1e10;
    dist = set_less(dist, distanceToRay(o, k, v1), nacb::Vec3d(1,0,0), vweights);
    dist = set_less(dist, distanceToRay(o, k, v2), nacb::Vec3d(0,1,0), vweights);
    dist = set_less(dist, distanceToRay(o, k, v3), nacb::Vec3d(0,0,1), vweights);

    nacb::Matrix M(3,2);
    nacb::Matrix b(3,1);
    for(int i=0; i<3; i++){
      M(i,0) = e1.data[i];
      M(i,1) = -k.data[i];
      
      b[i] = o.data[i] - v1.data[i];
    }
    nacb::Matrix x = nacb::Matrix::LlinLeastSq(M, b);
        
    if(x[0]>=0 && x[0]<=1){
      dist = set_less(dist, ((e1*x[0] + v1) - (o + k*x[1])).len(), nacb::Vec3d(1.0-x[0], x[0], 0.0), vweights);
    }

    for(int i=0; i<3; i++)M(i,0) = e2.data[i];
    
    x = nacb::Matrix::LlinLeastSq(M, b);
    
    if(x[0]>=0 && x[0]<=1){
      dist = set_less(dist, ((e2*x[0] + v1) - (o + k*x[1])).len(), nacb::Vec3d(1.0-x[0], 0.0, x[0]), vweights);
    }

    for(int i=0; i<3; i++){
      M(i,0) = e3.data[i];
      b[i] = o.data[i] - v2.data[i];
    }
    
    x = nacb::Matrix::LlinLeastSq(M, b);
    
    if(x[0]>=0 && x[0]<=1){
      dist = set_less(dist, ((e3*x[0] + v2) - (o + k*x[1])).len(), nacb::Vec3d(0.0, 1.0-x[0], x[0]), vweights);
    }
    return dist;
  }

  /*
    Return the minimum distance from this triangle to the point.
  */
  double distance(const nacb::Vec3f & point,
		  nacb::Vec3d * vweights = 0) const {
    double dist = 1e10;
    dist = set_less(dist, (v1 - point).len(), nacb::Vec3d(1,0,0), vweights);
    dist = set_less(dist, (v2 - point).len(), nacb::Vec3d(0,1,0), vweights);
    dist = set_less(dist, (v3 - point).len(), nacb::Vec3d(0,0,1), vweights);

    nacb::Vec3f proj = point - n*(point.dot(n) - d);
    nacb::Vec2f  p(i1.dot(proj), i2.dot(proj));
    p.x = std::max(p.x, 0.f);
    p.y = std::max(p.y, 0.f);
    
    p /= (p.x + p.y);
        
    nacb::Vec3d wt(p.x, p.y, 1.0-p.x-p.y);
    dist = set_less(dist, (v1 + e1*p.x + e2*p.y - point).len(), wt, vweights);

    return dist;
  }


  /*
    Intersect the ray with the triangle.  Return true if intersection happens.
    When intersected (and if not null) tout will be set to ray distance, and 
    a, b, c are set to bary centric weights (FIXME: should change these to be a vector).
  */
  bool intersect(const nacb::Vec3f & o,
		 const nacb::Vec3f & k,
		 double * tout = 0,
		 double * a = 0, 
		 double * b = 0, 
		 double * c = 0) const {
    
    double t = (d - o.dot(n))/(k.dot(n));

    nacb::Vec3f proj = (o + t*k) - v1;
    nacb::Vec2f p2(i1.dot(proj), i2.dot(proj));

    double sum = p2.x + p2.y;


    if(tout)*tout = t;

    const double eps = 1e-5;
    if(p2.x>=-eps && p2.y>=-eps && p2.x<=(1.0+eps) && p2.y<=(1.0+eps) && sum<=(1.0+eps)){
      // nacb::Vec3f rec = e1*p2.x + e2*p2.y + v1;
      // nacb::Vec3f diff = (o + t*k) - rec;

      if(a && b && c){
        static int log_message = 0;
        if (log_message < 10) {
          fprintf(stderr, "WARNING: %s:%d a,b,c were wrong before. Make sure code that uses them still works\n",
                  __FILE__, __LINE__);
          log_message++;
        }
	/*
	*a = p2.x/sum;
	*b = p2.y/sum;
	*c = (1.0 - *a - *b);
	*/
	*b = p2.x;
	*c = p2.y;
	*a = (1.0 - *b - *c);
      }
      return true;
    }
    return false;
  }
};

#endif // IK_TRIANGLE_H
