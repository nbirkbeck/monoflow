#include "Math.h"

extern "C" void dposv_(const char &uplo,               // (input)
		       const int &n,                   // (input)
		       const int &nrhs,                // (input)
		       double *a,                      // a[n][lda] (input/output)
		       const int &lda,                 // (input)
		       double *b,                      // b[nrhs][ldb] (input/output)
		       const int &ldb,                 // (input)
		       int &info                       // (output)
		       );




bool isNormal(nacb::Matrix & m){
  for(int i=0; i<m.n*m.m; i++)
    if(isinf(m.data[i]) || isnan(m.data[i]))
      return false;

  return true;
}


nacb::Mat5x5 unpack_sym_mat5x5(const nacb::Imagef & image, int x, int y){
  nacb::Mat5x5 m;

  int ind = 0;
  for(int i=0; i<5; i++){
    for(int j=i; j<5; j++){
      m(i, j) = image(x, y, ind);
      m(j, i) = m(i, j);
      ind++;
    }
  }
  return m;
}


void solve_symmetric_nxn(Matrix & M, Matrix & b, Matrix & x){
  if (M.m != M.n) 
    throw std::runtime_error("M must be square");

  if (x.m != M.m)
    x = Matrix(M.m, 1);

  int info = 0;
  memcpy(x.data, b.data, sizeof(double)*M.m);
  dposv_('U', M.m, 1, M.data,  M.m, x.data, M.m, info);
}


void solve_symmetric_4x4(Matrix & M, Matrix & b, Matrix & x){
  solve_symmetric_nxn(M, b, x);
}


Imagef get_D_phi_S(const nacb::Imagef U_grad[4], 
		   double hx, double hy,
                   double (*d_phi_S)(double)){
  nacb::Imagef D_phi_S(U_grad[0].w, U_grad[0].h, 2);

  for(int y=0; y<D_phi_S.h; y++){
    for(int x=0; x<D_phi_S.w; x++){ 
      D_phi_S(x, y, 0) = d_phi_S(sqr(U_grad[0](x, y, 0)/hx) + sqr(U_grad[0](x, y, 1)/hy));
      D_phi_S(x, y, 1) = d_phi_S(sqr(U_grad[1](x, y, 0)/hx) + sqr(U_grad[1](x, y, 1)/hy) +
				 sqr(U_grad[2](x, y, 0)/hx) + sqr(U_grad[2](x, y, 1)/hy) +
				 sqr(U_grad[3](x, y, 0)/hx) + sqr(U_grad[3](x, y, 1)/hy));
    }
  }
  return D_phi_S;
}
