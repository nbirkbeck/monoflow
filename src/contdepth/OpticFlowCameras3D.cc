/**
   Make it so you can use displace_offset_basis for optic flow.
*/
#include "autorecon/clbfile.h"
#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>
#include "ik/mesh.h"

int main(int ac, char * av[]) {
  if (ac < 2) {
    std::cerr << "Need the input sequence name." << std::endl;
    return -1;
  }

  std::string fileName = (boost::format(av[1]) % 0).str();
  std::string branch =  boost::filesystem::path(fileName).branch_path().string();
  std::string calib = branch + "/calib.clb";

  nacb::Image8 image;
  if (!image.read(fileName.c_str())) {
    std::cerr << "Cannot read image name." << std::endl;
    return -1;
  }

  double A[9] = {(double)image.w, 0, image.w/2.,
		 0, (double)image.w, image.h/2.,
		 0, 0, 1.};

  double E[16] = {1, 0, 0, 0,
		  0, 1, 0, 0,
		  0, 0, 1, 0,
		  0, 0, 0, 1};

  ClbFile::calibration_t clb(A, E);  
  clb.write(calib.c_str());

  double a = double(image.w)/double(image.h);
  Mesh mesh;
  mesh.vert.push_back(nacb::Vec3f(-0.5, -0.5/a, 1));
  mesh.vert.push_back(nacb::Vec3f( 0.5, -0.5/a, 1));
  mesh.vert.push_back(nacb::Vec3f( 0.5,  0.5/a, 1));
  mesh.vert.push_back(nacb::Vec3f(-0.5,  0.5/a, 1));

  mesh.tvert.push_back(nacb::Vec3f( 0, 1, 0));
  mesh.tvert.push_back(nacb::Vec3f( 1, 1, 0));
  mesh.tvert.push_back(nacb::Vec3f( 1, 0, 0));
  mesh.tvert.push_back(nacb::Vec3f( 0, 0, 0));

  mesh.restVert = mesh.vert;

  mesh.tris.push_back(Mesh::Tri());

  mesh.tris.back().vi[0] = 0;
  mesh.tris.back().vi[1] = 1;
  mesh.tris.back().vi[2] = 2;

  mesh.tris.back().tci[0] = 0;
  mesh.tris.back().tci[1] = 1;
  mesh.tris.back().tci[2] = 2;

  mesh.tris.push_back(Mesh::Tri());

  mesh.tris.back().vi[0] = 0;
  mesh.tris.back().vi[1] = 2;
  mesh.tris.back().vi[2] = 3;

  mesh.tris.back().tci[0] = 0;
  mesh.tris.back().tci[1] = 2;
  mesh.tris.back().tci[2] = 3;

  mesh.saveobj((branch + "/proxy.obj").c_str(), false, true);
  
		 
	    
		 

  return 0;
}
