/////////////////////////////////////////////////////////////////////////////////
//// 
////  Prototypes and definitions for sparse bundle adjustment
////  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
///////////////////////////////////////////////////////////////////////////////////

#ifndef _SBA_H_
#define _SBA_H_

#ifdef __cplusplus
extern "C" {
#endif

#define SBA_APPEND_UNDERSCORE_SUFFIX // undef this for AIX

#define SBA_MIN_DELTA     1E-06 // finite differentiation minimum delta
#define SBA_DELTA_SCALE   1E-04 // finite differentiation delta scale

#define SBA_OPTSSZ        4
#define SBA_INFOSZ        10
#define SBA_INIT_MU       1E-03
#define SBA_STOP_THRESH   1E-12
#define SBA_CG_NOPREC     0
#define SBA_CG_JACOBI     1
#define SBA_CG_SSOR       2
#define SBA_VERSION       "1.3 (Jun. 2006)"


/* Sparse matrix representation using Compressed Row Storage (CRS) format.
 * See http://www.netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000
 */

struct sba_crsm{
    int nr, nc;   /* #rows, #cols for the sparse matrix */
    int nnz;      /* number of nonzero array elements */
    int *val;     /* storage for nonzero array elements. size: nnz */
    int *colidx;  /* column indexes of nonzero elements. size: nnz */
    int *rowptr;  /* locations in val that start a row. size: nr+1.
                   * By convention, rowptr[nr]=nnz
                   */
};

/* sparse LM */

/* simple drivers */
extern int
sba_motstr_levmar(const int n, const int m, const int mcon, char *vmask, double *p, const int cnp, const int pnp, double *x, const int mnp,
           void (*proj)(int j, int i, double *aj, double *bi, double *xij, void *adata),
           void (*projac)(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata),
           void *adata, int itmax, int verbose, double opts[SBA_OPTSSZ], double info[SBA_INFOSZ]);

extern int
sba_mot_levmar(const int n, const int m, const int mcon, char *vmask, double *p, const int cnp, double *x, const int mnp,
           void (*proj)(int j, int i, double *aj, double *xij, void *adata),
           void (*projac)(int j, int i, double *aj, double *Aij, void *adata),
           void *adata, int itmax, int verbose, double opts[SBA_OPTSSZ], double info[SBA_INFOSZ]);

extern int
sba_str_levmar(const int n, const int m, char *vmask, double *p, const int pnp, double *x, const int mnp,
           void (*proj)(int j, int i, double *bi, double *xij, void *adata),
           void (*projac)(int j, int i, double *bi, double *Bij, void *adata),
           void *adata, int itmax, int verbose, double opts[SBA_OPTSSZ], double info[SBA_INFOSZ]);


/* expert drivers */
extern int
sba_motstr_levmar_x(const int n, const int m, const int mcon, char *vmask, double *p, int const cnp, const int pnp, double *x, const int mnp,
           void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
           void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
           void *adata, int itmax, int verbose, double opts[SBA_OPTSSZ], double info[SBA_INFOSZ]);

extern int
sba_mot_levmar_x(const int n, const int m, const int mcon, char *vmask, double *p, int const cnp, double *x, const int mnp,
           void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
           void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
           void *adata, int itmax, int verbose, double opts[SBA_OPTSSZ], double info[SBA_INFOSZ]);

extern int
sba_str_levmar_x(const int n, const int m, char *vmask, double *p, const int pnp, double *x, const int mnp,
           void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
           void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
           void *adata, int itmax, int verbose, double opts[SBA_OPTSSZ], double info[SBA_INFOSZ]);


/* interfaces to LAPACK routines: solution of linear systems, matrix inversion */
extern int sba_Axb_QR(double *A, double *B, double *x, int m, int iscolmaj);
extern int sba_Axb_QRnoQ(double *A, double *B, double *x, int m, int iscolmaj);
extern int sba_Axb_Chol(double *A, double *B, double *x, int m, int iscolmaj);
extern int sba_Axb_LU(double *A, double *B, double *x, int m, int iscolmaj);
extern int sba_Axb_SVD(double *A, double *B, double *x, int m, int iscolmaj);
extern int sba_Axb_BK(double *A, double *B, double *x, int m, int iscolmaj);
extern int sba_Axb_CG(double *A, double *B, double *x, int m, int niter, double eps, int prec, int iscolmaj);
extern int sba_mat_invert_LU(double *A, double *B, int m);
extern int sba_mat_invert_Chol(double *A, double *B, int m);

/* CRS sparse matrices manipulation routines */
extern void sba_crsm_alloc(struct sba_crsm *sm, int nr, int nc, int nnz);
extern void sba_crsm_free(struct sba_crsm *sm);
extern int sba_crsm_elmidx(struct sba_crsm *sm, int i, int j);
extern int sba_crsm_row_elmidxs(struct sba_crsm *sm, int i, int *vidxs, int *jidxs);
extern int sba_crsm_col_elmidxs(struct sba_crsm *sm, int j, int *vidxs, int *iidxs);

#ifdef __cplusplus
}
#endif

#endif /* _SBA_H_ */
