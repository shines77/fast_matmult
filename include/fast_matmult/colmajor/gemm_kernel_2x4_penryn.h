
#ifndef _GEMM_KERNEL_2X4_PENRYN_H_
#define _GEMM_KERNEL_2X4_PENRYN_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

void __CDECL _GCC_CDECL
gemm_kernel_2x4_penryn(const int m, const int n, const int k,
                       const cblas_float alpha,
                       const cblas_float *a, const int lda,
                       const cblas_float *b, const int ldb,
                       const cblas_float beta,
                       cblas_float *c, const int ldc,
                       const int offset);

#ifdef __cplusplus
}
#endif

#endif  /* _GEMM_KERNEL_2X4_PENRYN_H_ */
