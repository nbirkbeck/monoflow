/**
   FlowBasis3D.

   This file holds the base (and derived classes) to represent
   3D flow using a basis.  In fact, the classes represent a displacement
   and a 3D offset over time with a basis function. 

   Each basis type interprets data packed in an image as the coefficients
   of the basis.  The reason for this representation is so it is easy to
   upsample/downsample the coefficients for the multigrid implementation.
   
   TODO: 
   -Merge: getBasisMatrix, and derivative, can probably be merged somehow.
   -Create a way to instantiate the classes.

   (c) Neil Birkbeck, Jan 2011
*/
#ifndef FLOW_BASIS_3D
#define FLOW_BASIS_3D

#include <nmath/matrix.h>
#include <string>
#include <nimage/image.h>
#include <nmath/vec3.h>
#include <boost/filesystem/path.hpp>
#include <iostream>
#include <fstream>




/** \brief The base class used to represent the flow.
 */
class FlowBasis3D {
 public:
  typedef boost::shared_ptr<FlowBasis3D> ptr;

  virtual ~FlowBasis3D() { }

  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) = 0;

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) = 0;

  //!\brief Get the basis matrix for a set of times.
  virtual nacb::Matrix getBasisMatrix(std::vector<int>& times) = 0;

  //!\brief Get the number of elements in the basis.
  virtual int getNumDimensions() const  = 0;
  
  //!\brief Evaluate the displacement and offsets using the current coefficients.
  virtual nacb::Vec3f operator()(double & displace, float x, float y, const nacb::Imagef& D, int t) = 0;

  //!\brief Compute the derivative of the displacement/offset as a function of the current coefficients.
  virtual nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) = 0;

  virtual bool hasDisplacementComponent() const {
    return true;
  }

  //!\brief Some basis types (e.g., perFrame) are sparse and can save storage by packing differently.
  virtual bool isSparse() { return false; }

  //!\brief Return the interaction between variables at a given time (used to determine which coefficients must be stored)
  virtual nacb::Matrix getSparseInteractions(int t) {return nacb::Matrix(); }

  virtual std::string getType() const {
    return "standard";
  }


  /** \brief Extract a displacement/offset map from a 3D flow basis.
  */
  inline nacb::Imagef getDisplacementMap(const nacb::Imagef& D, int t) 
  {
    nacb::Imagef dispMap(D.w, D.h, 4);
    
    for (int y=0; y<dispMap.h; y++) {
      for (int x=0; x<dispMap.w; x++) {
	double displace;
	nacb::Vec3f offs = this->operator()(displace, x, y, D, t);
	
	dispMap(x, y, 0) = displace;
	dispMap(x, y, 1) = offs.x;
	dispMap(x, y, 2) = offs.y;
	dispMap(x, y, 3) = offs.z;
      }
    }
    
    return dispMap;
  }


  inline static FlowBasis3D::ptr create(const std::string& fileName, 
					nacb::Imagef& imageRepresentation);
};


class NullBasis3D: public FlowBasis3D {
 public:
  NullBasis3D() { }
  int getNumDimensions() const { return 4; }

  nacb::Matrix getBasisMatrix(std::vector<int>& times) {
    nacb::Matrix basisMatrix = nacb::Matrix::zeros(4, 0);
    return basisMatrix;
  }

  nacb::Vec3f operator()(double & displace, float x, float y, const nacb::Imagef& D, int t) {
    displace = 0;
    return nacb::Vec3f(0, 0, 0);
  }

  virtual nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    return nacb::Matrix();
  }

  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    return true;
  }

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    return true;
  }
};


/** \brief Standard constant velocity basis function.
 */
class StandardFlowBasis3D: public FlowBasis3D {
 public:
  StandardFlowBasis3D(int _referenceTime = 0): m_referenceTime(_referenceTime) { }

  int getNumDimensions() const { return 4; }

  nacb::Matrix getBasisMatrix(std::vector<int>& times) {
    nacb::Matrix basisMatrix = nacb::Matrix::zeros(4, 4*times.size());
    for (int j=0; j<(int)times.size(); j++) {
      basisMatrix(0, 4*j) = 1;
      basisMatrix(1, 4*j + 1) = (times[j] - m_referenceTime);
      basisMatrix(2, 4*j + 2) = (times[j] - m_referenceTime);
      basisMatrix(3, 4*j + 3) = (times[j] - m_referenceTime);
    }
    return basisMatrix;
  }

  nacb::Vec3f operator()(double & displace, float x, float y, const nacb::Imagef& D, int t) {
    displace = D.bilinear(x, y, 0);
    return nacb::Vec3f(D.bilinear(x, y, 1)*(t - m_referenceTime),
		       D.bilinear(x, y, 2)*(t - m_referenceTime),
		       D.bilinear(x, y, 3)*(t - m_referenceTime));
  }
  
  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Matrix  deriv = nacb::Matrix::eye(4, 4);
    deriv(1, 1) *= (t - m_referenceTime);
    deriv(2, 2) *= (t - m_referenceTime);
    deriv(3, 3) *= (t - m_referenceTime);
    return deriv;
  }

  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string leaf =  boost::filesystem::path(fileName).leaf().string();
    std::ofstream ofs(fileName.c_str());
    ofs << "basis3D standard 1 " << m_referenceTime << " " << (leaf + ".rfi");
    imageRepresentation.save((fileName + ".rfi").c_str());
    return true;
  }

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string branch =  boost::filesystem::path(fileName).branch_path().string();
    std::ifstream ifs(fileName.c_str());
    if (!ifs.good()) return false;
    std::string magic, type;
    int version;

    ifs >> magic >> type >> version;
    if (magic != "basis3D") return false;
    if (type != "standard") return false;
    if (version != 1) return false;

    std::string localFile;
    m_referenceTime = 0;
    ifs >> m_referenceTime >> localFile;

    imageRepresentation.load((branch + "/" + localFile).c_str());

    if (imageRepresentation.nchannels != 4) 
      return false;
    return true;
  }

  virtual std::string getType() const {
    return "standard";
  }

 protected:
  int m_referenceTime;
};



/** \brief Standard constant velocity basis function.
 */
class StandardFlowBasis2D: public FlowBasis3D {
 public:
  StandardFlowBasis2D(int _referenceTime = 0): m_referenceTime(_referenceTime) { }

  int getNumDimensions() const { return 2; }

  nacb::Matrix getBasisMatrix(std::vector<int>& times) {
    nacb::Matrix basisMatrix = nacb::Matrix::zeros(2, 4*times.size());
    for (int j=0; j<(int)times.size(); j++) {
      basisMatrix(0, 4*j + 1) = (times[j] - m_referenceTime);
      basisMatrix(1, 4*j + 2) = (times[j] - m_referenceTime);
    }
    return basisMatrix;
  }

  nacb::Vec3f operator()(double & displace, float x, float y, const nacb::Imagef& D, int t) {
    displace = 0;
    return nacb::Vec3f(D.bilinear(x, y, 0)*(t - m_referenceTime),
		       D.bilinear(x, y, 1)*(t - m_referenceTime),
		       0);
  }
  
  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Matrix  deriv = nacb::Matrix::zeros(4, 2);
    deriv(1, 0) = (t - m_referenceTime);
    deriv(2, 1) = (t - m_referenceTime);
    return deriv;
  }

  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string leaf =  boost::filesystem::path(fileName).leaf().string();
    std::ofstream ofs(fileName.c_str());
    ofs << "basis3D standard2D 1 " << m_referenceTime << " " << (leaf + ".rfi");
    imageRepresentation.save((fileName + ".rfi").c_str());
    return true;
  }

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string branch =  boost::filesystem::path(fileName).branch_path().string();
    std::ifstream ifs(fileName.c_str());
    if (!ifs.good()) return false;
    std::string magic, type;
    int version;

    ifs >> magic >> type >> version;
    if (magic != "basis3D") return false;
    if (type != "standard2D") return false;
    if (version != 1) return false;

    std::string localFile;
    m_referenceTime = 0;
    ifs >> m_referenceTime >> localFile;

    imageRepresentation.load((branch + "/" + localFile).c_str());

    if (imageRepresentation.nchannels != 2) 
      return false;
    return true;
  }

  virtual bool hasDisplacementComponent() const {
    return false;
  }

  virtual std::string getType() const {
    return "standard2D";
  }

 protected:
  int m_referenceTime;
};



/** \brief Per frame basis 3D.

    The reference frame doesn't have parameters for offset, so there
    is some index adjustements in this file.
    
    The basis for 3 frames (where 0 is reference) would be:
    basis = [1, 0, 0, 0,  1, 0, 0, 0,  1, 0, 0, 0;
             0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0;
	     0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0;
	     0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0;
	     0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0;
	     0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 1, 0;
	     0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 1];
	     
 */
class PerFrameFlowBasis3D: public FlowBasis3D {
 public:
  PerFrameFlowBasis3D(int _referenceTime = 0, int _numTime = 4):
   m_referenceTime(_referenceTime), m_numTime(_numTime) { }
  
  int getNumDimensions() const {return (m_numTime - 1)*3 + 1; }

  nacb::Vec3f operator()(double & displace, 
			 float x, float y, 
			 const nacb::Imagef& D,
			 int t) {
    displace = D((int)x, (int)y, 0);
    if (1 + 3*(t - 1) >= D.nchannels || t < 0) {
      return nacb::Vec3f(0, 0, 0);
    }
    if (t == m_referenceTime)
      return nacb::Vec3f(0, 0, 0);

    int ti = (t >= m_referenceTime) ? (t - 1) : t;

    return nacb::Vec3f(D((int)x, (int)y, 1 + 3*ti),
		       D((int)x, (int)y, 2 + 3*ti),
		       D((int)x, (int)y, 3 + 3*ti));
  }

  virtual nacb::Matrix getBasisMatrix(std::vector<int>& times) {
    nacb::Matrix basis = nacb::Matrix::zeros((m_numTime - 1)*3 + 1, 4*times.size());
    
    for (int i=0; i<(int)times.size(); ++i)
      basis(0, 4*i) = 1;

    for (int i=0; i<min(m_numTime, (int)times.size()); i++) {
      int ri = i;
      if (i == m_referenceTime) continue;
      if (i > m_referenceTime)
	ri--;

      basis(1 + 3*ri, 4*i + 1) = 1;
      basis(2 + 3*ri, 4*i + 2) = 1;
      basis(3 + 3*ri, 4*i + 3) = 1;
    }
    return basis;
  }

  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Matrix deriv = nacb::Matrix::zeros(4, (m_numTime - 1)*3 + 1);

    deriv(0, 0) = 1;
    if (t == m_referenceTime) return deriv;

    if (t < m_numTime) {
      int ti = t;
      if (t >= m_referenceTime) 
	ti = t - 1;

      deriv(1, 1 + 3*ti) = 1;
      deriv(2, 2 + 3*ti) = 1;
      deriv(3, 3 + 3*ti) = 1;
    }
    return deriv;
  }

  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string leaf =  boost::filesystem::path(fileName).leaf().string();
    std::ofstream ofs(fileName.c_str());
    ofs << "basis3D perframe 1 " << m_referenceTime << " " << m_numTime << " " << (leaf + ".rfi");

    imageRepresentation.save((fileName + ".rfi").c_str());
    return true;
  }

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string branch =  boost::filesystem::path(fileName).branch_path().string();
    std::ifstream ifs(fileName.c_str());
    if (!ifs.good()) return false;
    std::string magic, type;
    int version;
    
    ifs >> magic >> type >> version;

    if (magic != "basis3D") return false;
    if (type != "perframe") return false;
    if (version != 1) return false;

    std::string localFile;
    ifs >> m_referenceTime >> m_numTime >> localFile;
    
    std::cout << "ref:" << m_referenceTime << " " <<  m_numTime << " " << localFile << std::endl;
    imageRepresentation.load((branch + "/" + localFile).c_str());

    if (imageRepresentation.nchannels != (m_numTime - 1)*3 + 1) 
      return false;
    return true;
  }

  virtual bool isSparse() { return true; }

  virtual nacb::Matrix getSparseInteractions(int t) {
    nacb::Matrix m = nacb::Matrix::zeros(getNumDimensions(), 1);
    m[0] = 1;

    int ti = t;
    if (t == m_referenceTime) 
      return m;

    if (t >= m_referenceTime) 
      ti = t - 1;
    
    if (3*ti + 3 < getNumDimensions()) {
      m[3*ti + 1] = 1;
      m[3*ti + 2] = 1;
      m[3*ti + 3] = 1;
    }
    return m;
  }

  virtual std::string getType() const {
    return "perframe";
  }

 protected:
  int m_referenceTime;
  int m_numTime;
};


class CosineBasis {
 public:
  CosineBasis(int _numTime = 0,
	      int _referenceTime = 0): 
  m_numTime(_numTime), 
    m_referenceTime(_referenceTime) { }

  double evaluateBasis(int basisIndex, int t) {
    return cos(M_PI/m_numTime * (0.5 + t - m_referenceTime)*basisIndex);
  }

 protected:
  int m_numTime;
  int m_referenceTime;
};


/** \brief Cosine basis for offset in z component only (assumes it is aligned with normal so it
           can also model displacement).
    
 */
class CosineDepthOffset: public FlowBasis3D, protected CosineBasis {
 public:
 CosineDepthOffset(int _referenceTime = 0, int _numBasis = 4, int _numTime = 4): 
  CosineBasis(_numTime, _referenceTime), m_numBasis(_numBasis) { }

  int getNumDimensions() const {return m_numBasis; }

  nacb::Vec3f operator()(double & displace, 
			 float x, float y, 
			 const nacb::Imagef& D,
			 int t) {
    nacb::Vec3f dx(0, 0, 0);
    displace = 0;

    for (int i=0; i<D.nchannels; i++) {
      double value = getWeight(t) *  evaluateBasis(i, t) * D((int)round(x), (int)round(y), i);
      dx.data[2] += value;
    }
    return dx;
  }

  virtual nacb::Matrix getBasisMatrix(std::vector<int>& times) {
    nacb::Matrix basis = nacb::Matrix::zeros(getNumDimensions(), 4*times.size());
    for (int t=0; t<(int)times.size(); t++) {
      for (int i=0; i<m_numBasis; i++) {
	int basisIndex = i;
	basis(i, 4*t + 3) = getWeight(times[t]) * evaluateBasis(basisIndex, t);
      }
    }
    return basis;
  }

  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Matrix deriv = nacb::Matrix::zeros(4, getNumDimensions());

    for (int i=0; i<D.nchannels; i++) {
      double value = getWeight(t) * evaluateBasis(i, t);
      deriv(3, i) = value;
    }
    return deriv;
  }
  
  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string leaf =  boost::filesystem::path(fileName).leaf().string();
    std::ofstream ofs(fileName.c_str());
    ofs << "basis3D cosinedepth 1 " << m_referenceTime << " " 
	<< m_numTime << " " << m_numBasis << " " << (leaf + ".rfi");;

    imageRepresentation.save((fileName + ".rfi").c_str());
    return true;
  }

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string branch =  boost::filesystem::path(fileName).branch_path().string();
    std::ifstream ifs(fileName.c_str());
    if (!ifs.good()) return false;
    std::string magic, type;
    int version;
    
    ifs >> magic >> type >> version;

    if (magic != "basis3D") return false;
    if (type != "cosinedepth") return false;
    if (version != 1) return false;

    std::string localFile;
    ifs >> m_referenceTime >> m_numTime >> m_numBasis >> localFile;

    imageRepresentation.load((branch + "/" + localFile).c_str());

    if (imageRepresentation.nchannels != getNumDimensions()) {
      std::cout << "channels not equal:" << imageRepresentation.nchannels << " " << getNumDimensions() << std::endl;
      return false;
    }
    return true;
  }

  //!\brief Use this weight to zero the offset at the reference time.
  double getWeight(int t) const {
    return 1.0;
  }

  //!\breif This displaces along the offset, so it doesn't have a displacement component.
  virtual bool hasDisplacementComponent() const {
    return false;
  }


 protected:
  int m_numBasis;
};



/** \brief Cosine basis 3D.
 */
class CosineFlowBasis3D: public FlowBasis3D, protected CosineBasis {
 public:
 CosineFlowBasis3D(int _referenceTime = 0, int _numBasis = 4, int _numTime = 4, bool _useDepthBasis = false): 
  CosineBasis(_numTime, _referenceTime), m_numBasis(_numBasis), m_depthBasis(_useDepthBasis) { }
  
  int getNumDimensions() const {return m_numBasis * 3 + m_depthBasis; }

  nacb::Vec3f operator()(double & displace, 
			 float x, float y, 
			 const nacb::Imagef& D,
			 int t) {
    nacb::Vec3f dx(0, 0, 0);

    if (m_depthBasis)
      displace = D((int)round(x), (int)round(y), 0);
    else
      displace = 0;

    for (int i=m_depthBasis; i<D.nchannels; i++) {
      int basisIndex = int((i - m_depthBasis)/3);
      double value = getWeight(t, (i - m_depthBasis) % 3) *  evaluateBasis(basisIndex, t) * D((int)round(x), (int)round(y), i);
      dx.data[(i - m_depthBasis) % 3] += value;
    }
    return dx;
  }

  virtual nacb::Matrix getBasisMatrix(std::vector<int>& times) {
    nacb::Matrix basis = nacb::Matrix::zeros(3*m_numBasis + m_depthBasis, 4*times.size());
    for (int t=0; t<(int)times.size(); t++) {
      if (m_depthBasis)
	basis(0, 4*t) = 1.0;

      for (int i=0; i<m_numBasis*3; i++) {
	int basisIndex = int(i/3);
	basis(i + m_depthBasis, 4*t + (i % 3) + 1) = getWeight(times[t], i % 3) * evaluateBasis(basisIndex, t);
      }
    }
    return basis;
  }

  nacb::Matrix derivative(float x, float y, const nacb::Imagef& D, int t) {
    nacb::Matrix deriv = nacb::Matrix::zeros(4, getNumDimensions());

    if (m_depthBasis)
      deriv(0, 0) = 1;

    for (int i=m_depthBasis; i<D.nchannels; i++) {
      int basisIndex = int((i - m_depthBasis)/3);
      double value = getWeight(t, (i - m_depthBasis)%3) * evaluateBasis(basisIndex, t);
      deriv(((i - m_depthBasis) % 3) + 1, i) = value;
    }
    return deriv;
  }

  virtual bool save(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string leaf =  boost::filesystem::path(fileName).leaf().string();
    std::ofstream ofs(fileName.c_str());
    ofs << "basis3D cosine 2 " << m_referenceTime << " " 
	<< m_numTime << " " << m_numBasis << " " << (leaf + ".rfi") << " " << m_depthBasis;

    imageRepresentation.save((fileName + ".rfi").c_str());
    return true;
  }

  virtual bool load(const std::string& fileName,
		    nacb::Imagef& imageRepresentation) {
    std::string branch =  boost::filesystem::path(fileName).branch_path().string();
    std::ifstream ifs(fileName.c_str());
    if (!ifs.good()) return false;
    std::string magic, type;
    int version;
    
    ifs >> magic >> type >> version;

    if (magic != "basis3D") return false;
    if (type != "cosine") return false;
    if (version > 2) return false;

    // The default was to have a depth basis.
    m_depthBasis = true;

    std::string localFile;
    ifs >> m_referenceTime >> m_numTime >> m_numBasis >> localFile;
    
    if (version >= 2) 
      ifs >> m_depthBasis;

    std::cout << m_referenceTime << " " << m_numTime << " " << m_numBasis << " " << localFile << std::endl;
    imageRepresentation.load((branch + "/" + localFile).c_str());

    if (imageRepresentation.nchannels != getNumDimensions()) {
      std::cout << "channels not equal:" << imageRepresentation.nchannels << " " << getNumDimensions() << std::endl;
      return false;
    }
    return true;
  }

  virtual std::string getType() const {
    return "cosine";
  }

  //!\brief Use this weight to zero the offset at the reference time.
  double getWeight(int t, int co) const {
    if (m_depthBasis || co != 2)
      return  t != m_referenceTime;
    else  // The z coordinate always has an element.
      return 1.0;
  }

  virtual bool hasDisplacementComponent() const {
    return m_depthBasis;
  }

 protected:
  int m_numBasis;
  bool m_depthBasis; //!< Whether we have a separate displacement component
};


inline FlowBasis3D::ptr FlowBasis3D::create(const std::string& fileName, 
					    nacb::Imagef& imageRepresentation) 
{
  
  {
    FlowBasis3D::ptr basis(new StandardFlowBasis3D);
    if (basis->load(fileName, imageRepresentation))
      return basis;
  }

  {
    FlowBasis3D::ptr basis(new CosineFlowBasis3D);
    if (basis->load(fileName, imageRepresentation))
      return basis;
  }

  {
    FlowBasis3D::ptr basis(new PerFrameFlowBasis3D);
    if (basis->load(fileName, imageRepresentation))
      return basis;
  }  

  {
    FlowBasis3D::ptr basis(new CosineDepthOffset);
    if (basis->load(fileName, imageRepresentation))
      return basis;
  }  

  {
    FlowBasis3D::ptr basis(new StandardFlowBasis2D);
    if (basis->load(fileName, imageRepresentation))
      return basis;
  }  
  return FlowBasis3D::ptr();
}




#endif // FLOW_BASIS_3D
