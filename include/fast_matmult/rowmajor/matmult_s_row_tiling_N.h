
#ifndef _MATMULT_S_ROW_TILING_N_H_
#define _MATMULT_S_ROW_TILING_N_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <stdlib.h>

#include <fast_matmult/fast_matmult.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_N_sse2_2x4.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_N_sse2_2x4_packed.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_N_sse2_4x1.h>
#include <fast_matmult/colmajor/gemm_kernel_2x4_penryn.h>

#ifdef __cplusplus
extern "C" {
#endif

void matmult_s_row_tiling_MxK_N(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

void matmult_s_row_tiling_KxM_N(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

void matmult_s_row_tiling_k_n_m_KxM_N(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_N_H_ */
