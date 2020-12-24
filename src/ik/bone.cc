#pragma implementation "bone.h"
#include "bone.h"

double Bone::cons_scales[5]={1,1,1,1,1};

bool Bone::drawBounds=false;
int Bone::updated=0;

/*
bool operator==(const Bone::bone_affects_t & b1,
		const Bone::bone_affects_t & b2){
  return b1.first==b2.first && b1.second==b2.second;
}
*/
