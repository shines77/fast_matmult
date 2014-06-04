
#ifndef _MATMULT_S_ROW_TILING_K_SSE2_2x4_PACKED_H_
#define _MATMULT_S_ROW_TILING_K_SSE2_2x4_PACKED_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

void matmult_s_row_tiling_NxM_K_sse2_2x4_packed(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_K_SSE2_2x4_PACKED_H_ */
