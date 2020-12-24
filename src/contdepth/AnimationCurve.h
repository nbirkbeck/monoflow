#ifndef ANIMATION_CURVE_H
#define ANIMATION_CURVE_H


#include <map>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>

class AnimationCurve {
public:

  AnimationCurve() {
    name = "";
  }

  AnimationCurve(const std::string & name) {
    this->name = name;
  }

  const double & get(double t) const {
    std::map<double, double>::const_iterator low = keys.lower_bound(t);

    if(low != keys.end())
      return low->second;

    static double res = 0;
    return res;
  }

  double & get(double t) {
    std::map<double, double>::iterator low = keys.lower_bound(t);

    if(low != keys.end())
      return low->second;

    static double res = 0;
    return res;
  }

  const double & operator()(double t) const {
    return get(t);
  }

  const double & operator[](double t) const {
    return get(t);
  }

  double & operator()(double t){
    return get(t);
  }

  double & operator[](double t){
    return get(t);
  }

  void insert(double t, double value) {
    keys[t] = value;
  }

  double max() const {
    if (!keys.size()) return 0.0;
    return keys.rbegin()->first;
  }

  std::map<double, double> keys;
  std::string name;

 private:
  template<class Archive>
    void save(Archive & ar, unsigned int version) const
    {
      
      ar & (const std::string &)name;
      ar & (const size_t &)(keys.size());

      std::map<double, double>::const_iterator it = keys.begin();
      std::map<double, double>::const_iterator last = keys.end();

      printf("Saving %s\n", name.c_str());

      while (it != last) {
	ar & (const double &) it->first;
	ar & (const double &) it->second;
	it++;
      }
    }

  template<class Archive>
    void load(Archive & ar, const unsigned int version)
  {
    ar & name;

    size_t size;
    ar & size;

    for (int i=0; i<(int)size; i++) {
      double key, value;
      ar & key & value;
      
      keys[key] = value;
    }
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER();
  friend class boost::serialization::access;
};


class AnimationCurves {
public:
  AnimationCurves() {

  }

  bool count(const std::string & name) {
    return curves.count(name);
  }

  AnimationCurve & operator[](const std::string & name) {
    if (!curves.count(name)){
      printf("no curve.\n");
      curves[name] = AnimationCurve(name);
    }
    return curves[name];
  }

  void print(std::ostream & o) const {
    o << "count: " << curves.size();
  }

  double max() const {
    std::map<std::string, AnimationCurve>::const_iterator it(curves.begin());
    std::map<std::string, AnimationCurve>::const_iterator last(curves.end());

    double val = 0.0;
    for (; it != last; it++) {
      val = std::max(val, it->second.max());
    }   
    return val;
  }
  
  std::map<std::string, AnimationCurve> curves;

 private:
  template<class Archive>
    void save(Archive & ar, const unsigned int version) const 
    {
      size_t sz = curves.size();
      ar & sz;
      
      std::map<std::string, AnimationCurve>::const_iterator it(curves.begin());
      
      for(; it != curves.end(); it++){
	printf("outputing curves.\n");
	const AnimationCurve & curve = it->second;
	ar & (const AnimationCurve &) curve;
      }
    }

  template<class Archive>
    void load(Archive & ar, const unsigned int version)
  {
    size_t sz;
    
    ar & sz;

    for (int i=0; i<(int)sz; i++) {
      AnimationCurve curve;

      ar & curve;
      
      curves[curve.name] = curve;
    }    
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER();
  friend class boost::serialization::access;
};





namespace std {
  
  inline std::ostream & operator<<(std::ostream & o, const AnimationCurves & curves) {
    curves.print(o);
    return o;
  }
  
};


#endif
