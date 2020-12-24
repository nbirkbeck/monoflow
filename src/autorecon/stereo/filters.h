#ifndef STEREO_FILTERS_H
#define STEREO_FILTERS_H

#include <nimage/image.h>

nacb::Imagef fast_median_reject(const nacb::Imagef & im,
				int hs, 
				int x1, int x2,
				float reject);

nacb::Imagef median_reject(const nacb::Imagef & im, 
			   int hs = 7, float reject = 0);
nacb::Imagef trilateral(const nacb::Imagef & im, const nacb::Imagef & sc,
			int hs=7, float reject=0);

void fill_holes(nacb::Imagef & im, float reject, int maxsize);

#endif //STEREO_FILTERS_H
