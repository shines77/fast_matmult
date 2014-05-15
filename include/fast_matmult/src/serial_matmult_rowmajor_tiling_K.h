
#ifndef _SERIAL_MATMULT_ROWMAJOR_TILING_K_H_
#define _SERIAL_MATMULT_ROWMAJOR_TILING_K_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/fast_matmult.h>

#include <stdlib.h>

#if (_WIN32 | _WIN64) && _MSC_VER
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

namespace annlab {

void serial_matmult_rowmajor_tiling_K_NxM(
    const int M, const int N, const int K,
    const float_t alpha,
    const float_t *A, const int lda,
    const float_t *B, const int ldb,
    const float_t beta,
    float_t *C, const int ldc);

}  /* namespace annlab */

#endif  /* _SERIAL_MATMULT_ROWMAJOR_TILING_K_H_ */
