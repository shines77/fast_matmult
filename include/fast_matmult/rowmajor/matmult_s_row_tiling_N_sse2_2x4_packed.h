
#ifndef _MATMULT_S_ROW_TILING_N_SSE2_2X4_PACKED_H_
#define _MATMULT_S_ROW_TILING_N_SSE2_2X4_PACKED_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

void matmult_s_row_tiling_KxM_N_sse2_2x4_packed(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_N_SSE2_2X4_PACKED_H_ */
