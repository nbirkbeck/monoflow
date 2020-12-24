#ifndef OPTIC_FLOW_H
#define OPTIC_FLOW_H

#include <nimage/image.h>
#include <nmath/matrix.h>

/*! \brief warp_back image im1 using the displacements in disp.
*/
nacb::Imagef flow_warp_back(const nacb::Imagef & disp, const nacb::Imagef & im1);


/*! \brief Compute the fundamental matrix from the flow field, using the mask
           to determine which pixels should be used.  p1out and p2out will contain
	   the point sets used (when non-null).
*/
nacb::Matrix flow_fundamental_matrix(const nacb::Imagef & disp, const nacb::Image8 & mask,
				     nacb::Matrix * p1out = 0, nacb::Matrix * p2out = 0);


/*! \brief Return the length of the flow field (as an image).
*/
nacb::Imagef flow_length(const nacb::Imagef & disp);




/*! \brief compose two flow-fields (linear update, from grid-based paper, \sa flow_grid.cc).
   Returns the compositional combination of two displacement fields U and update (upd)

   Notice that the compositional methods updates the flow by transforming the update 
   by looking up the displacement from existing displacement using newly displaced
   coordinates and then adding to this the updated displacement.
*/
nacb::Imagef flow_compose(const nacb::Imagef & flo, const nacb::Imagef & update);

bool flow_save_length(const nacb::Imagef & disp, const char * name);

nacb::Imagef flow_read_flo(const char * fname);
bool flow_save_flo(const nacb::Imagef & image, const char * fname);


bool flow_save_crossfade(const char * basename,
			 const nacb::Imagef & im1,
			 const nacb::Imagef & im2,
			 const nacb::Imagef & d1,
			 const nacb::Imagef & d2,
			 int nt = 10);

std::vector<nacb::Imagef> flow_crossfade(const nacb::Imagef & im1,
					 const nacb::Imagef & im2,
					 const nacb::Imagef & d1,
					 const nacb::Imagef & d2,
					 int nt = 10);

nacb::Imagef opticFlow(const nacb::Imagef & im1,const nacb::Imagef & im2,int hwin,
		       int x0, int x1, int y0, int y1);
nacb::Imagef opticFlowTraditional(const nacb::Imagef & im1, const nacb::Imagef & im2, const nacb::Imagef & estimate, int hwin, nacb::Imagef * det = 0);

#endif //OPTIC_FLOW_H
