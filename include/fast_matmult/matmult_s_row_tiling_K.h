
#ifndef _MATMULT_S_ROW_TILING_K_H_
#define _MATMULT_S_ROW_TILING_K_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <stdlib.h>

#if (_WIN32 | _WIN64) && _MSC_VER
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

#include <fast_matmult/fast_matmult.h>

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

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_TILING_K_H_ */
