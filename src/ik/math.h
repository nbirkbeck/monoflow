#ifndef IK_MATH_H
#define IK_MATH_H

#include <boost/function.hpp>
#include "armature.h"

class OptimizerLogger {
 public:
  enum Level {
    OptimizationLevel = 0,
    IterationLevel = 1,
    SubIterationLevel = 2,
  };

  OptimizerLogger(){ }

  bool iterationComplete() const {
    return true;
  }

  void operator()(const Level & level,
		  const std::string & message) const {
    return log(level, message);
  }

  void operator()(const std::string & message, 
		  const Level & level = OptimizationLevel) const {
    return log(level, message);
  }

  virtual void log(const Level & level, 
		   const std::string & message) const { 
    std::cout << message;
  }

  virtual void warning(const std::string & message) const { 
    std::cerr << message;
  }

  virtual void error(const std::string & message) const { 
    std::cerr << message;
  }
};

class NullLogger : public OptimizerLogger {
 public:
  NullLogger(){

  }
  virtual void log(const Level & level, 
		   const std::string & message) const { }
};


class Basis1D {
 public:
  std::vector<double (*)(double )> funcs;

  Basis1D(){
    funcs.push_back(cb0);
    funcs.push_back(cb1);
    funcs.push_back(cb2);
    funcs.push_back(cb3);
  }

  double eval(double t, double coeff[]){
    double val = 0;
    for(int i=0; i<funcs.size(); i++)
      val += funcs[i](t)*coeff[i];
    return val;
  }

  double basis_eval(int fi, double t){
    return funcs[fi](t);
  }
  
  static double cb0(double t){
    return t*t*t;
  }

  static double cb1(double t){
    return 3*t*t*(1.0-t);
  }

  static double cb2(double t){
    return 3*t*(1.0-t)*(1.0-t);
  }

  static double cb3(double t){
    return (1.0-t)*(1.0-t)*(1.0-t);
  }

  static double d_cb0(double t){
    return 3*t*t;
  }

  static double d_cb1(double t){
    return 3*t*t*(1.0-t);
  }

  static double d_cb2(double t){
    return 3*t*(1.0-t)*(1.0-t);
  }

  static double d_cb3(double t){
    return (1.0-t)*(1.0-t)*(1.0-t);
  }
};

inline double uniform_random(double n, double x){
  return drand48()*(x-n)+n;
}


//FIXME: in libs/um_lbfgs/lbfgs.o
extern "C" void lbfgs_(const int & n,
		       const int & m,
		       double * X,
		       double & F,
		       double * G,
		       const int & diagco,
		       double * diag,
		       int * iprint,
		       const double & eps,
		       const double & xtol,
		       double * W,
		       int & iflag);

//FIXME: this is defined in vsphereset.cc
double golden_section(boost::function<double (double)> func,
		      double t,double e);

double xorFuncEvalAsym(const Image8 & psil, const Image8 & sils);
double xorFuncEval8(const Image8 & psil,const Image8 & sils);

double jaccard(const Image8& psil, const Image8 & sils, int numViews);

Mat4x4 getProjection4x4(Mat3x3 & K, Mat4x4 & E);

Matrix powell(const boost::function<double (const Matrix &)> & fnc, Matrix & t, int maxits = 10, 
	      const OptimizerLogger & logger = OptimizerLogger());

void powell(const boost::function<double (Armature *,Mesh *,Matrix &,Matrix & ,double)> & fnc,
	    Armature * arm,Mesh * mesh);

void graddesc(const boost::function<double (Armature * ,Mesh *,Matrix &,Matrix & ,double)> & fnc,Armature * arm,Mesh * mesh);

void lbfgs(const boost::function<double (Armature *,Mesh *,Matrix &,Matrix & ,double)> & fnc,Armature * arm,Mesh * mesh, double fdspc=0.0125);


//Allocates and returns the gradient.
Matrix fd_gradient(const boost::function<double (Matrix & )> & func, Matrix & x, double dx);

//Does not allocate unless necessary.
double func_and_fd_gradient(const boost::function<double (Matrix & )> & func, Matrix & x, double dx, Matrix & gout);

void lbfgs(const boost::function<double (Matrix & parms, Matrix & gradient)> &, Matrix & x);

#endif
