#ifndef FLOW_VARIATIONAL_H
#define FLOW_VARIATIONAL_H

#include <nimage/image.h>

struct robust_stats {
  double moved, res, len, erel;
  robust_stats(){  }

  robust_stats(double _moved, double _res, double _len, double _erel){
    moved = _moved;
    res = _res;
    len = _len;
    erel = _erel;
  }

  void print(const char * header = "", int pads = 0){
    for(int i=0; i<pads; i++)putchar(' ');
    printf("%s:  moved=%lf, residual=%lf, length=%lf, rel=%lf\n", header, moved, res, len, erel);
  }
};


double phi_1en6(double s2);
double phi_1en4(double s2);
double phi_1en2(double s2);

double d_phi_const(double s2);
double d_phi_1en6(double s2);
double d_phi_1en4(double s2);
double d_phi_1en2(double s2);

nacb::Imagef get_D_phi_S(const nacb::Imagef & D_grad, double hx, double hy, double (*d_phi_S)(double));

#endif // FLOW_VARIATIONAL_H
