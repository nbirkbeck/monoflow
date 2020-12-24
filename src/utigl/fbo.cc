#include "fbo.h"

int checkFrameBufferStatus(){
  GLenum status;
  status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch(status) {						  
  case GL_FRAMEBUFFER_COMPLETE_EXT:			
    break;							  
  case GL_FRAMEBUFFER_UNSUPPORTED_EXT:				  
    /* choose different formats */				  
    printf("framebuffer unsupported\n");
    break;							  
  default:
    printf("framebuffer error\n");
    printf("exiting\n");
    exit(0);
  }
  return 0;
}
