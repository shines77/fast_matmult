
#ifndef _MATMULT_S_ROW_TILING_K_H_
#define _MATMULT_S_ROW_TILING_K_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_K_sse2_2x4.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_K_sse2_2x4_packed.h>

#ifdef __cplusplus
extern "C" {
#endif

void matmult_s_row_tiling_NxM_K(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

// ========================================================================

void matmult_s_row_tiling_MxN_K_transB(unsigned int M, unsigned int K, unsigned int N,
                                       float_t *A, float_t *B, float_t *C);

// ========================================================================

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_K_H_ */
