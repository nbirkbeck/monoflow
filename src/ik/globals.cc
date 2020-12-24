#include "globals.h"
#include <GL/gl.h>

using namespace nacb;

nacb::Imagef unprojectDepth(nacb::Imagef & z, double minz, double maxz, bool origintl){
  nacb::Imagef depth(z.w, z.h, 1);

  for(int y=0; y<z.h; y++){
    for(int x=0; x<z.w; x++){
      double zn = z(x, y)*2.0 - 1.0;
      //during projection:
      //zn = (zw * mat[10] + mat[14])/(zw * mat[11]);
      //zw  = mat[14]/(zn * (mat[11] -  mat[10]));

      double mat10 = (origintl?-1: 1)*(minz+maxz)/(minz-maxz);
      double mat11 = (origintl? 1:-1);      
      double mat14 =2*maxz*minz/(minz-maxz);

      depth(x, y) = mat14/(zn * mat11 - mat10);
    }
  }
  return depth;
}


void myProjection(double fx,double fy,
                  double px,double py,
                  int imw,int imh,
                  double minz,double maxz,
		  int origintl)
{
  double mat[16];
  memset(mat,0,sizeof(mat));

  mat[0]=fx*2/imw;
  
  //Aspect
  //mat[4]=intrinsics[4]*2/imw;
  mat[5]=fy*2/imh;

  mat[8]=px*2/imw-1;
  mat[9]=py*2/imh-1;
  mat[10]=(origintl?-1: 1)*(minz+maxz)/(minz-maxz);
  mat[11]=(origintl? 1:-1);
  
  mat[14]=2*maxz*minz/(minz-maxz);

  glMultMatrixd(mat);
}

bool tag_eq(void * t1, const char * t2){
  return *((uint32_t*)t1)==*((uint32_t*)t2);
}

bool readCLB(const char * fname,
	     Mat3x3 & K,
	     Mat4x4 & E){
  FILE * file=fopen(fname,"rb");
  if(!file)return false;
    
  uint32_t tag_len[2];
  if(2!=fread(&tag_len,4,2,file)){
    fclose(file);
    return false;
  }
  if(!tag_eq(tag_len,"CLB ")){
    fclose(file);
    return false;
  }
  int got=0;
  while(2==fread(&tag_len,4,2,file)){
    double m[16];
    if(tag_eq(tag_len,"INTR")){
      size_t r = fread(m,sizeof(double),9,file);
      if (r != 9) {
        fprintf(stderr, "Error reading INTR\n");
      }
      for(int i=0;i<9;i++)
	K(i/3,i%3)=m[i];
      got++;
    }
    else if(tag_eq(tag_len,"EXTR")){
      size_t r = fread(m,sizeof(double),16,file);
      if (r != 16) {
        fprintf(stderr, "Error reading EXTR\n");
      }
      for(int i=0;i<16;i++)
	E(i/4,i%4)=m[i];
      got++;
    }
    else
      fseek(file,tag_len[1],SEEK_CUR);
  }
  fclose(file);
  return got==2;
}


void drawArrow(const Vec4d & from,const Vec4d & to,double aspect,double head){
  glPushMatrix();

  double mod[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,mod);

  Vec3d yd,zd,xd=Vec3d(to.x,to.y,to.z)-Vec3d(from.x,from.y,from.z);
  double len=xd.normalize();
  
  zd=Vec3d(mod[2],mod[6],mod[10]);
  zd.normalize();
  yd=zd.cross(xd);
  yd.normalize();
  zd=xd.cross(yd);

  mod[0]=xd.x;
  mod[1]=xd.y;
  mod[2]=xd.z;
  mod[3]=0;

  mod[4]=yd.x;
  mod[5]=yd.y;
  mod[6]=yd.z;
  mod[7]=0;

  mod[8]= zd.x;
  mod[9]= zd.y;
  mod[10]=zd.z;
  mod[11]=0;

  mod[12]=from.x;
  mod[13]=from.y;
  mod[14]=from.z;
  mod[15]=1;
  glMultMatrixd(mod);

  glScalef(len,len*aspect,len*aspect);
  
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(head,0,0);
  glVertex3f(0,-0.1,0);
  glVertex3f(head,-0.1,0);
  glVertex3f(head,-0.2,0);
  glVertex3f(1.0,0.0,0);
  glVertex3f(head, 0.2,0);
  glVertex3f(head, 0.1,0);
  glVertex3f(0, 0.1,0);
  glVertex3f(0,-0.1,0);
  
  glEnd();
  glPopMatrix();
}
