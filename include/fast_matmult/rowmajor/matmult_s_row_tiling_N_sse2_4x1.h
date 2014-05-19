
#ifndef _MATMULT_S_ROW_TILING_N_SSE2_4X1_H_
#define _MATMULT_S_ROW_TILING_N_SSE2_4X1_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

void __CDECL
matmult_s_row_tiling_N_sse2_4x1(const int m, const int n, const int k,
                                const float_t *alpha,
                                const float_t *a, const int lda,
                                const float_t *b, const int ldb,
                                const float_t *beta,
                                float_t *c, const int ldc);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_N_SSE2_4X1_H_ */
