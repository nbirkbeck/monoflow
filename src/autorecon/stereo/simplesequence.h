#ifndef SIMPLE_SEQUENCE_H
#define SIMPLE_SEQUENCE_H

#include <nimage/image.h>
#include <nmath/matrix.h>
#include "../clbfile.h"
#include <boost/format.hpp>
#include <libgen.h>

template <class T>
class BaseSequence {
 public:
  nacb::Image<T>   image;
  nacb::Matrix A, E;
  nacb::Matrix P;
  nacb::Matrix kc;
  nacb::Vec3d cc;

  BaseSequence(){ 
    A = nacb::Matrix::eye(3, 3);
    E = nacb::Matrix::eye(4, 4);
    P = nacb::Matrix::eye(3, 4);
    cc = nacb::Vec3d(0, 0, 0);
    kc = nacb::Matrix::zeros(5, 1);
  }

  BaseSequence(const char * filename){ 
    A = nacb::Matrix::eye(3, 3);
    E = nacb::Matrix::eye(4, 4);
    P = nacb::Matrix::eye(3, 4);
    kc = nacb::Matrix::zeros(5, 1);

    if (!initCalibrationFromFile(filename)) {
      A = nacb::Matrix::eye(3, 3);
      E = nacb::Matrix::eye(4, 4);
      P = nacb::Matrix::eye(3, 4);
      kc = nacb::Matrix::zeros(5, 1);

      updateDerivedCameraInfo();
    }
  }

  virtual ~BaseSequence() { }
  
  virtual bool initCalibrationFromFile(const char * filename) {
    ClbFile::flags_t fags = ClbFile::read(filename,
					  A.data, E.data, P.data, kc.data);
    bool ret =((fags & ClbFile::EXTRINSICS) && (fags & ClbFile::INTRINSICS));
    if (ret)
      updateDerivedCameraInfo();
    // printf("Internals for camera:%d\n", i);
    // A.printMatlab("A");
    // E.printMatlab("E");
    // d.printMatlab("d");
    // printf("\n");
    return ret;
  }

  virtual void updateDerivedCameraInfo() {
    nacb::Matrix Einv = E.inverse();
    cc = nacb::Vec3d(Einv(0,3), Einv(1,3), Einv(2,3));
    P = A*nacb::Matrix::eye(3,4)*E;
  }  

  virtual bool load(int index) {
    return false;
  }
};


template <class T>
class ImageSequence : public BaseSequence<T> {
public:

  using BaseSequence<T>::image;
  using BaseSequence<T>::A;
  using BaseSequence<T>::P;
  using BaseSequence<T>::E;
  using BaseSequence<T>::cc;
  
  std::string basename;
  bool staticCamera;

  static double s_subfrac;
  static int s_subdiv;

  ImageSequence(): BaseSequence<T>(), staticCamera(true) { }

  ImageSequence(const std::string & basename, int i=0){
    this->basename = basename;    
    staticCamera = true;

    char name[1024];
    snprintf(name, 1024, basename.c_str());

    const char * dir = dirname(name);
    std::string calibName = std::string(dir + std::string("/calib.clb"));
    
    if (!BaseSequence<T>::initCalibrationFromFile(calibName.c_str())) {
      staticCamera = false;
      
      calibName = ((boost::format(basename) % i).str() + ".clb");
      if (!BaseSequence<T>::initCalibrationFromFile(calibName.c_str()))
	throw std::string("cannot dynamic calibration file: ") + calibName;
    }

    if(!load(0)){
      throw std::string("cannot load first image");
    }

    updateDerivedCameraInfo();
  }
  
  virtual void updateDerivedCameraInfo() {
    nacb::Matrix Einv = E.inverse();

    for(int i=0; i<s_subdiv; i++){
      nacb::Matrix sc = nacb::Matrix::eye(3,3);
      sc(0,0) = sc(1,1) = s_subfrac;
      A = sc*A;
    }

    cc = nacb::Vec3d(Einv(0,3), Einv(1,3), Einv(2,3));
    P = A*nacb::Matrix::eye(3,4)*E;
  }

  bool load(int index){
    char fname[1024];
    snprintf(fname, 1024, basename.c_str(), index);
    
    image = nacb::Image<T>(fname);
    
    for(int i=0; i<s_subdiv; i++){
      image = image.resize((int)(s_subfrac*image.w), (int)(s_subfrac*image.h));
    }

    // Non-static cameras need to reload themselves
    if (!staticCamera) {
      std::string fileName = ((boost::format(basename) % index).str() + ".clb");
      if (!BaseSequence<T>::initCalibrationFromFile(fileName.c_str()))
	return false;

      std::cerr << "Warning: non-static camera has not been tested. If it works, remove this message: " << __FILE__ << "," << __LINE__ << "\n";
      updateDerivedCameraInfo();
    }

    return image.w*image.h*image.nchannels>1;
  }

  void draw(){
    
  }
};

#define SimpleSequence ImageSequence


#endif // SIMPLE_SEQUENCE_H
