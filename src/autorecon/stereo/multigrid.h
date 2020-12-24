/** \brief Helper code for multigrid methods.
 */
#ifndef MULTIGRID_H
#define MULTIGRID_H

#include <nimage/image.h>
#include <nmath/matrix.h>
#include <assert.h>
#include <vector>
#include <nmath/vec2.h>

namespace mg {


  /** \brief perform a restriction/prolongation using image resize operations.
   */
  struct restriction_image {
    int w, h, neww, newh;
  
    
    restriction_image(int _w = 0, int _h = 0, 
		      int _neww = 0, int _newh = 0) : w(_w), h(_h), neww(_neww), newh(_newh){
		    }
  
    nacb::Imagef apply(const nacb::Imagef & image){
      return image.resize(neww, newh);
    }

    nacb::Imagef prolong(const nacb::Imagef & image){
      return image.resize(w, h);
    }
  };



  /** \brief the averaging restriction/prolongation used by brox (I think).
      Somewhat slower than restriction_image.
  */
  struct restriction_average {

    //!\brief Yet another way to store a sparse matrix.
    struct entry {
      int      x, y;
      double    val;
    
    entry(int _x=0, int _y=0,
	  double _val=0.0) : 
      x(_x), y(_y), val(_val){ }

    entry(const nacb::Vec2<int> & _loc,
	  double _val) : x(_loc.x), y(_loc.y), val(_val){ }
    };

    std::vector<entry> * locoeff;
    std::vector<entry> * hicoeff;
    int w, h, neww, newh;

    restriction_average(): locoeff(0), hicoeff(0), w(0), h(0), neww(0), newh(0) {}
    restriction_average(int _w, int _h, 
			int _neww, int _newh);

    ~restriction_average();

    void construct(int _w, int _h, int _neww, int _newh);

    nacb::Imagef prolong(const nacb::Imagef & input);
    nacb::Imagef apply(const nacb::Imagef & input);

    void print();
  private:
    void operator=(const restriction_average & other);
  };


  struct wcycle_solver {
    int  presmooth_its;
    int postsmooth_its;
    int     ngrids;
    int    wcycles;
    bool nonlinear;
    int   base_its;

    wcycle_solver(int _ngrids = 3,
		  int _wcycles = 3,
		  int _pre = 10, 
		  int _post = 10,
		  bool _nonlinear  = false) : 
      presmooth_its(_pre), postsmooth_its(_post),
      ngrids(_ngrids), wcycles(_wcycles), nonlinear(_nonlinear) 
    { 	
      base_its = presmooth_its + postsmooth_its;
    }

    template <typename problem_t>
    void wcycle_iteration(problem_t & problem, int depth = 0){
      if(depth<=0 || !problem.can_restrict()){
	for(int i=0; i<base_its; i++)
	  problem.smooth();
	return;
      }

      for(int cyclei=0; cyclei<wcycles; cyclei++){
	for(int i=0; i<presmooth_its; i++)
	  problem.smooth();
      
	problem_t restricted = problem.restriction();

	auto H_orig_sol = restricted.copy_solution();

	wcycle_iteration(restricted, depth - 1);
		
	if(!problem.is_linear())
	  problem.add_to_solution(problem.prolong(restricted.get_solution() - H_orig_sol));
	else
	  problem.add_to_solution(problem.prolong(restricted.get_solution()));
	
	for(int i=0; i<postsmooth_its; i++)
	  problem.smooth();
      }
    }

    template<class problem_t>
    void solve(problem_t & problem, int topits = 1){
      for(int its = 0; its<topits; its++){
	wcycle_iteration(problem, ngrids);
      }
    }
  };

  struct basic_solver {
    template<class problem_t>
    void solve(problem_t & problem, int topits = 10000){
      for(int its = 0; its<topits; its++){
	problem.smooth();
      }
    }
  };
};


#endif //MULTIGRID_H
