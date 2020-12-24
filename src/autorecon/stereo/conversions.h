#ifndef STEREO_CONVERSIONS_H
#define STEREO_CONVERSIONS_H

#include <nimage/image.h>
#include <nmath/matrix.h>



nacb::Imagef disparityToDepth(const nacb::Imagef & disp,
			      const nacb::Matrix & Pn1,
			      const nacb::Matrix & Pn2,
			      double minx1, double minx2,
			      double disp_reject, double depth_reject);

void writeRawDepth(const char * imname, const nacb::Imagef & im);
nacb::Imagef readRawDepth(const char * imname);
#endif
