
#ifndef _MATMULT_S_ROW_TILING_K_SSE2_2x4_H_
#define _MATMULT_S_ROW_TILING_K_SSE2_2x4_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

void matmult_s_row_tiling_NxM_K_sse2_2x4(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

void matmult_s_row_tiling_MxN_K_transB_sse2_2x4(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_K_SSE2_2x4_H_ */
