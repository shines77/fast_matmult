
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

void __CDECL
gemm_kernel_2x4_penryn(const int m, const int n, const int k,
                       const float_t alpha,
                       const float_t *a, const int lda,
                       const float_t *b, const int ldb,
                       const float_t beta,
                       float_t *c, const int ldc,
                       const int offset);

#ifdef __cplusplus
}
#endif

#endif  /* _GEMM_KERNEL_2X4_PENRYN_H_ */
