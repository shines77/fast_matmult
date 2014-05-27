
#ifndef _MATMULT_BLAS_DGEMM_H_
#define _MATMULT_BLAS_DGEMM_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/cblas_def.h>

#ifdef __cplusplus
extern "C" {
#endif

void cblas_dgemm(const Cblas_Order order, const Cblas_Transpose transa,
                 const Cblas_Transpose transb,
                 const cblas_int m, const cblas_int n, const cblas_int k,
                 const float_t *alpha,
                 const float_t *a, const cblas_int lda,
                 const float_t *b, const cblas_int ldb,
                 const float_t *beta,
                 float_t *c, const cblas_int ldc,
                 cblas_func gemm_func);

#ifdef __cplusplus
}
#endif

#endif /* _MATMULT_BLAS_DGEMM_H_ */
