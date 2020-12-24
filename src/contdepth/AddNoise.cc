#include <nmisc/commandline.h>
#include <nimage/image.h>

int main(int ac, char * av[]) {
  nacb::CommandLine cline;
  double noise = 1.0;

  cline.registerOption("noise", "The amount of noise to add", &noise, 'n');
  cline.parse(ac, av);

  for (int i=optind; i<ac; i++) {
    nacb::Image8 image;
    image.read(av[i]);
    
    for (int y=0; y<image.h; y++) {
      for (int x=0; x<image.w; x++) {
	for (int k=0; k<3; k++) {
	  float value = image(x, y, k);

	  value = std::max(0.f, std::min(255.0f, (float)(value + (2.0*drand48() - 1.0)*noise)));
	  image(x, y, k) = (unsigned char)value;
	}
      }
    }
    image.save(av[i]);
  }

  return 0;
}
