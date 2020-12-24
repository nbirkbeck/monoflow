#ifndef FBO_H
#define FBO_H

#include <GL/gl.h>
#include <GL/glut.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#ifndef glGenFramebuffers 
#define glGenFramebuffers glGenFramebuffersEXT
#endif

#ifndef glGenRenderbuffers
#define glGenRenderbuffers glGenRenderbuffersEXT
#endif

#ifndef glBindFramebuffer
#define glBindFramebuffer glBindFramebufferEXT
#endif

#ifndef glBindRenderbuffer
#define glBindRenderbuffer glBindRenderbufferEXT
#endif

#ifndef glFramebufferTexture2D
#define glFramebufferTexture2D glFramebufferTexture2DEXT
#endif

#ifndef glRenderbufferStorage
#define glRenderbufferStorage glRenderbufferStorageEXT
#endif

#ifndef  glFramebufferRenderbuffer
#define glFramebufferRenderbuffer glFramebufferRenderbufferEXT
#endif 

#ifndef  glDeleteRenderbuffers
#define glDeleteRenderbuffers glDeleteRenderbuffersEXT
#endif 

#ifndef  glDeleteFramebuffers
#define glDeleteFramebuffers glDeleteFramebuffersEXT
#endif 

int checkFrameBufferStatus();

class FrameBufferObject{
 protected:

  GLuint fb;
  GLuint depth;
  GLuint stencil;
  GLuint color;
  int width,height;
 public:
  enum{
    USE_COLOR=0x1,
    USE_DEPTH=0x2,
    USE_STENCIL=0x4
  };
  //no option for turning on or off a particular buffer
  FrameBufferObject(int w=256,int h=256,
		    int use_buffers=USE_COLOR|USE_DEPTH){
    fb=depth=stencil=color=0;
    width=height=0;
    glGenFramebuffers(1,&fb);
    if(use_buffers&USE_STENCIL){
      glGenRenderbuffers(1,&stencil);
      assert(stencil);
    }
    if(use_buffers&USE_DEPTH){
      glGenRenderbuffers(1,&depth);
      assert(depth);
    }
    if(use_buffers&USE_COLOR){
      glGenTextures(1,&color);
      assert(color);
    }
    assert(fb);
    resize(w,h);
  }
  int po2(int n){
    int val=16;
    while(val<n)val*=2;
    return val;
  }
  void resize(int neww,int newh){
    neww=po2(neww);
    newh=po2(newh);
    
    if(neww==width&&newh==height)return;
    printf("initializing fbo %d %d\n",neww,newh);

    glBindFramebuffer(GL_FRAMEBUFFER_EXT, fb);

    if(color){
      glBindTexture(GL_TEXTURE_2D,color);
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T, GL_CLAMP);
    
      unsigned char * data = new unsigned char[neww*newh*4];
      for(int i=0;i<neww*newh*4;i++)
	data[i]=rand()%256;
      glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,neww,newh,0,GL_RGBA,GL_UNSIGNED_BYTE,data);
      delete [] data;
      
      glFramebufferTexture2D(GL_FRAMEBUFFER_EXT,
			     GL_COLOR_ATTACHMENT0_EXT,
			     GL_TEXTURE_2D,color,0);
    }
    
    if(depth){
      glBindRenderbuffer(GL_RENDERBUFFER_EXT, depth);
      glRenderbufferStorage(GL_RENDERBUFFER_EXT,
			    GL_DEPTH_COMPONENT,neww,newh);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER_EXT,
				GL_DEPTH_ATTACHMENT_EXT,
				GL_RENDERBUFFER_EXT, depth);
    }
    if(stencil){
      printf("initing stencil, hasn't worked as of yet\n");
      glBindRenderbuffer(GL_RENDERBUFFER_EXT, stencil);
      glRenderbufferStorage(GL_RENDERBUFFER_EXT,
			    GL_STENCIL_INDEX,neww,newh);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER_EXT,
				GL_STENCIL_ATTACHMENT_EXT,
				GL_RENDERBUFFER_EXT,stencil);
    }
    width=neww;
    height=newh;
    glBindFramebuffer(GL_FRAMEBUFFER_EXT, 0);
  }
  ~FrameBufferObject(){
    if(color)glDeleteTextures(1,&color);
    if(depth)glDeleteRenderbuffers(1,&depth);
    if(stencil)glDeleteRenderbuffers(1,&stencil);
    if(fb)glDeleteFramebuffers(1,&fb);
    color=depth=stencil=fb=0;
  }
  void bind(int bindit){
    glBindFramebuffer(GL_FRAMEBUFFER_EXT,bindit?fb:0);
    checkFrameBufferStatus();
  }
  void useColorTexture(GLuint id){
    if(id==0){
      glFramebufferTexture2D(GL_FRAMEBUFFER_EXT,
			     GL_COLOR_ATTACHMENT0_EXT,
			     GL_TEXTURE_2D,color,0);
    }
    else{
      glFramebufferTexture2D(GL_FRAMEBUFFER_EXT,
			     GL_COLOR_ATTACHMENT0_EXT,
			     GL_TEXTURE_2D,id,0);
    }
  }
  int getWidth(){
    return width;
  }
  int getHeight(){
    return height;
  }
  GLuint getColorTexture(){
    return color;
  }
  GLuint getDepthTexture(){
    return depth;
  }
 private:
  //this is undefined because texture deletion
  void operator=(const FrameBufferObject & other){
  }
};


#endif
