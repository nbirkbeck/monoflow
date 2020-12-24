#ifndef OFFSET_MESH_H
#define OFFSET_MESH_H

#include <GL/gl.h>
#include <nimage/image.h>
#include <nmath/matrix.h>
#include "autorecon/stereo/depthmesh.h"


class OffsetMesh : public DepthMesh {
 public:
  OffsetMesh(const nacb::Image8 & image,
	     const nacb::Imagef & depth,
	     const nacb::Imagef & offset,
	     const nacb::Matrix & P);

  void draw();

  void setOffsetImage(const nacb::Imagef & offset);

  void setOffsetScale(float f){
    offsetScale = f;
  }

 protected:
  GLuint offsetTex;
  float offsetScale;

 private:
  void operator=(OffsetMesh & offsetMesh){ }
};


#endif // OFFSET_MESH_H
