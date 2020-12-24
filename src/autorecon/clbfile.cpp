/* 
   \file clbfile.cpp
   \brief Implementation file for reading in calibration information.  

   Copyright (c) 2009, Neil Birkbeck      
*/

#include "clbfile.h"
#include <stdint.h>
#include <stdlib.h>

namespace ClbFile{
  std::ostream & operator<<(std::ostream & o, const ClbFile::flags_t & f){
    o << "(";
    if(f==NOTHING)
      o << "NOTHING";
    if(f & INTRINSICS)
      o << "INTRINSICS ";
    if(f & EXTRINSICS)
      o << "EXTRINSICS ";
    if(f & PROJECTION)
      o << "PROJECTION ";
    if(f & DISTORTION)
      o << "DISTORTION ";
    return o << ")";
  }

  void print_matrix(std::ostream & o,const char * nm,
		    const double * data,int m,int n){
    o << nm << "= [...\n";
    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
	o << data[i*n+j] << " ";
      }
      if(i==m-1)	
	o << "];";
      o << std::endl;
    }
  }


  flags_t read(const char * fileName,
	       double * A, double * E,
	       double * P, double * d){
    FILE * file;

    struct tag_len_t{
      uint32_t tag;
      uint32_t len;
    }tag_len;

    flags_t res = NOTHING;


    file = fopen(fileName, "rb");

    if(!file)
      return res;
    
    if(!fread(&tag_len, sizeof(tag_len), 1, file)){
      fclose(file);
      return res;
    }

    if(ClbFile::tag_equal(tag_len.tag,"CLB ")){

      while( fread(&tag_len, sizeof(tag_len), 1, file) ){
	if( tag_equal(tag_len.tag, "INTR") && tag_len.len==9*8){
	  size_t r = fread(A, tag_len.len, 1, file);
          if (r != 1) {
            fprintf(stderr, "Error reading INTR\n");
          }
	  res = (flags_t) (res | INTRINSICS);
	}
	else if( tag_equal(tag_len.tag, "EXTR") && tag_len.len==16*8){
	  size_t r = fread(E, tag_len.len, 1, file);
          if (r != 1) {
            fprintf(stderr, "Error reading EXTR\n");
          }
	  res = (flags_t) (res | EXTRINSICS);
	}
	else if( tag_equal(tag_len.tag, "DIST") && tag_len.len==5*8){
	  size_t r = fread(d, tag_len.len, 1, file);
          if (r != 1) {
            fprintf(stderr, "Error reading DIST\n");
          }
	  res = (flags_t) (res | DISTORTION);
	}
	else if( tag_equal(tag_len.tag, "PROJ") && tag_len.len==12*8){
	  size_t r = fread(P, tag_len.len, 1, file);
          if (r != 1) {
            fprintf(stderr, "Error reading PROJ\n");
          }
	  res = (flags_t) (res | PROJECTION);
	}
	else { 
	  //Unknown tag, skip the bytes
	  fseek(file, tag_len.len, SEEK_CUR);
	}
      }
    }
    else {
      fclose(file);
      return res;
    }
    
    fclose(file);
    return res;
  }

  flags_t read(const char * fileName, calibration_t & t){
    t.flags = read(fileName, t.A, t.E, t.P, t.d);
    return t.flags;
  }
};

std::ostream & operator<<(std::ostream & o, const ClbFile::calibration_t & c){
  c.print(o);
  return o;
}
