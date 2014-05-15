
#ifndef _SERIAL_MATMULT_ROWMAJOR_TILING_N_SSE2_4X1_H_
#define _SERIAL_MATMULT_ROWMAJOR_TILING_N_SSE2_4X1_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "fast_matmult/common_asm.h"
#include <fast_matmult/fast_matmult.h>

namespace annlab {

void __CDECL
serial_matmult_rowmajor_tiling_N_sse2_4x1(const int m, const int n, const int k,
                                          const float_t *alpha,
                                          const float_t *a, const int lda,
                                          const float_t *b, const int ldb,
                                          const float_t *beta,
                                          float_t *c, const int ldc);

}  /* namespace annlab */

#endif  /* _SERIAL_MATMULT_ROWMAJOR_TILING_N_SSE2_4X1_H_ */
