/* 
   \file   clbfile.h
   \brief  Header for routines to read in calibration file.

   Copyright (c) 2009, Neil Birkbeck      
*/
#ifndef CLBFILE_H
#define CLBFILE_H

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

namespace ClbFile{
  enum flags_t{
    NOTHING=0x0,
    INTRINSICS=0x1,
    EXTRINSICS=0x2,
    PROJECTION=0x4,
    DISTORTION=0x8
  };

  std::ostream & operator<<(std::ostream & o, const ClbFile::flags_t & f);
  void print_matrix(std::ostream & o,const char * nm,
		    const double * data,int m,int n);

  struct calibration_t {
    double A[9]; //3x3 row major
    double E[16];//4x4 row major
    double P[12];//3x4 row major, should be identity
    double d[5]; //5-vector
    flags_t flags; //what this structure actually contains.

    calibration_t(){
      flags = NOTHING;
    }

    calibration_t(const double * intr,
		  const double * extr=0,
		  const double * dist=0){
      flags = (flags_t)INTRINSICS;
      if(extr)
	flags = (flags_t)(flags | EXTRINSICS);
      if(dist)
	flags = (flags_t)(flags | DISTORTION);
      
      memcpy(A, intr, sizeof(double)*9);
      
      if(extr)
	memcpy(E, extr, sizeof(double)*16);

      if(dist)
	memcpy(d, dist, sizeof(double)*5);
    }
    void print(std::ostream & o) const {
      o << flags << "\n";
      if(flags & INTRINSICS)
	print_matrix(o,"A",A,3,3);
      
      if(flags & EXTRINSICS)
	print_matrix(o,"E",E,4,4);
      
      if(flags & PROJECTION)
	print_matrix(o,"P",P,3,4);

      if(flags & DISTORTION)
	print_matrix(o,"d",d,1,5);
    };
    void write(const char * fname) const {
      int len = 
	((flags & INTRINSICS)? 8+9*8: 0) +
	((flags & EXTRINSICS)? 8+16*8: 0) +
	((flags & PROJECTION)? 8+12*8: 0) +
	((flags & DISTORTION)? 8+5*8: 0);

      FILE * file = fopen(fname, "wb");
      fwrite("CLB ", 4, 1, file);
      fwrite(&len, 4, 1, file);
      
      if(flags & INTRINSICS){
	len = 9*8;
	fwrite("INTR",4,1,file);
	fwrite(&len,4,1,file);
	fwrite(A,sizeof(double),9,file);
      }
      if(flags & EXTRINSICS){
	len = 16*8;
	fwrite("EXTR",4,1,file);
	fwrite(&len,4,1,file);
	fwrite(E,sizeof(double),16,file);
      }
      if(flags & PROJECTION){
	len = 12*8;
	fwrite("PROJ",4,1,file);
	fwrite(&len,4,1,file);
	fwrite(P,sizeof(double),12,file);
      }
      if(flags & DISTORTION){
	len = 5*8;
	fwrite("DIST",4,1,file);
	fwrite(&len,4,1,file);
	fwrite(d,sizeof(double),5,file);
      }
      fclose(file);
    }
  };

  inline bool tag_equal(const int& i1, const char * i2){
    return i1 == *((int *)i2);
  }
  
  /** \brief calibration from fileName.

      \param A  pointer to  9 doubles (row major intrinsics 3x3)
      \param E  pointer to 16 doubles (row major intrinsics 4x4)
      \param P  pointer to 12 doubles (row major intrinsics 3x4)
      \param d  pointer to 5 doubles (distortion coefficients
              same as camera calibration toolbox for matlab)

      The projection matrix is always eye(3,4), only included in case
      there was orthographic projection.  You can safely ignore.

      Distortion coefficients may all be zero.
      
      \return returns a set of flags denoting what was contained in the file.
      
  */
  flags_t read(const char * fileName,
	       double * A, double * E,
	       double * P, double * d);

  flags_t read(const char * fileName, calibration_t & t);
};


#endif
