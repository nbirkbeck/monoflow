#ifndef DISPARITY_SEQUENCE_H
#define DISPARITY_SEQUENCE_H

#include <nmath/matrix.h>
#include "autorecon/stereo/simplesequence.h"

class DisparitySequence: public SimpleSequence<unsigned char> {

public:
  DisparitySequence(){ }
  DisparitySequence(const char * str): SimpleSequence<unsigned char>(str) {

  }
  
  bool loadDisparity(int i){
    
    char fname[1024];

    char * dname = strdup(basename.c_str());
    snprintf(fname, 1024, "%s/disp-%04d.rfi", dirname(dname), i);
    disp.load(fname);
    free(dname);

    dname = strdup(basename.c_str());
    snprintf(fname, 1024, "%s/flow-%04d.rfi", dirname(dname), i);
    flow.load(fname);
    free(dname);

    return (disp.w == image.w && disp.h == image.h);
  }
  
  DisparitySequence getDownsampled(int downsample) const {
    DisparitySequence seq;

    seq.disp = disp.copy();
    seq.flow = flow.copy();
    seq.E = E.copy();
    seq.cc = cc;
    seq.basename = basename;
    seq.image = image.copy();
    seq.A = A.copy();

    nacb::Matrix ds = nacb::Matrix::eye(3, 3);
    ds(0, 0) = ds(1, 1) = 0.5;

    for(int i=0; i<downsample; i++){
      seq.A = ds * seq.A;
      seq.P = seq.A*nacb::Matrix::eye(3, 4)*E;
      const nacb::Image8 & img = seq.image;
      
      seq.disp = disp.resize(img.w/2, img.h/2);
      seq.flow = flow.resize(img.w/2, img.h/2);
      
      seq.image = img.resize(img.w/2, img.h/2);
    }

    return seq;
  }

  static std::vector<DisparitySequence> loadSequences(char * av[], int num, int downsample = 0){
    std::vector<DisparitySequence> seqs;
    
    try{
      for (int i=0; i<num; i++){
	seqs.push_back(DisparitySequence(av[i]));
	if(seqs.back().loadDisparity(0)){
	  printf("Has disparity %s\n", av[i]);
	}      
      }
    }
    catch(const std::string & str){
      std::cout << "Exception: " << str << "\n";
    }
    
    return seqs;
  }

  static std::vector<DisparitySequence> downsample(const std::vector<DisparitySequence> & seqs, 
						   int downsample = 0){
    std::vector<DisparitySequence> results;
    
    for(int si=0; si<(int)seqs.size(); si++){
      results.push_back(seqs[si].getDownsampled(downsample));   
    }
    
    return results;
  }

  nacb::Imagef disp;
  nacb::Imagef flow;
};



typedef DisparitySequence SeqType;

#endif
