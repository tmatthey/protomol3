#include <protomol/type/Real.h>
#ifndef PROTOMOL_LAPACK_H
#define PROTOMOL_LAPACK_H

namespace ProtoMol {
  namespace Lapack {
    bool isEnabled();
    void dgemv(char *transA, int *m, int *n, Real *alpha, Real *A,
               int *lda, Real *x, int *incx, Real *beta, Real *Y,
               int *incY);
    void dsyev(char *jobz, char *uplo, int *n, Real *a, int *lda,
               Real *w, Real *work, int *lwork, int *info);
    void dsyevr(char *jobz, char *range, char *uplo, int *n, Real *a,
                int *lda, Real *vl, Real *vu, int *il, int *iu,
                Real *abstol, int *m, Real *w, Real *z,  int *ldz,
                int *isuppz, Real *work, int *lwork, int *iwork,
                int *liwork, int *info);
    Real dlamch(char *cmach);
    void dgemm(char *transA, char *transB, int *m, int *n, int *k,
               Real *alpha, Real *A, int *lda, Real *B, int *ldb,
               Real *beta, Real *C, int *l);
    Real ddot(int *n, Real *x, int *incx, Real *y, int *incy);
    Real dnrm2(int *n, Real *x, int *incx);
    void dpotri(char *transA, int *n, Real *A, int *lda, int *info);
    void dpotrf(char *transA, int *n, Real *A, int *lda, int *info);
    void dposv(char *transA, int *n, int *nrhs, Real *a, int *lda,
               Real *b, int *ldb,int *info);
    void dtrmm(char *sideA, char *ulA, char *transA, char *diagA,
               int *m, int *n, Real *alpha, Real *A, int *lda,
               Real *B, int *ldb);
    void dtrsm(char *sideA, char *ulA, char *transA, char *diagA,
               int *m, int *n, Real *alpha, Real *A, int *lda,
               Real *B, int *ldb);
  }
}

#endif // PROTOMOL_LAPACK_H
