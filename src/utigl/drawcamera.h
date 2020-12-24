/**
   I got sick of writing this all the time.
   Needed to put this into a common place.

   Need: 
   -draw the camera frustum.
   -optionally draw an image.
   -selection.
   -scale.

   -camera
   drawCamera(P, );
*/

#ifndef UTIGL_DRAWCAMERA_H
#define UTIGL_DRAWCAMERA_H

#include <nmath/matrix.h>

using namespace nacb;

class DrawCamera {
 public:
  Matrix Kinv4x4;
  Matrix Einv;
  int w, h;

  double scale;
  bool enableTexture;

  DrawCamera(const Matrix & K, const Matrix & E, int _w, int _h) {
    scale = 1.0;
    w = _w;
    h = _h;
    enableTexture = false;

    Einv = E.inverse();

    Kinv4x4 = Matrix::eye(4, 4);
    Kinv4x4.setSubmatrix(0, 0, K.inverse());
  }

  void setEnableTexture(bool val) {
    enableTexture = val;
  }
  
  void draw() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glMultMatrixd((Einv * Kinv4x4).transpose().data);

    glBegin(GL_LINE_LOOP);
    glVertex3f(0, 0, scale);
    glVertex3f(w*scale, 0, scale);
    glVertex3f(w*scale, h*scale, scale);
    glVertex3f(0, h*scale, scale);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, scale);

    glVertex3f(0, 0, 0);
    glVertex3f(w*scale, 0, scale);

    glVertex3f(0, 0, 0);
    glVertex3f(w*scale, h*scale, scale);

    glVertex3f(0, 0, 0);
    glVertex3f(0, h*scale, scale);
    glEnd();

    if (enableTexture) {
      glEnable(GL_TEXTURE_2D);

      glBegin(GL_QUADS);
       glTexCoord2f(0, 0);
       glVertex3f(0, 0, scale);
       glTexCoord2f(1, 0);
       glVertex3f(w*scale, 0, scale);
       glTexCoord2f(1, 1);
       glVertex3f(w*scale, h*scale, scale);
       glTexCoord2f(0, 1);
       glVertex3f(0, h*scale, scale);
      glEnd();

      glDisable(GL_TEXTURE_2D);
    }

    glPopMatrix();
  }
};


#endif // UTIGL_DRAWCAMERA_H
