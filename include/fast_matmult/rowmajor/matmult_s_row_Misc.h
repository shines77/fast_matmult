
#ifndef _MATMULT_S_ROW_MISC_H_
#define _MATMULT_S_ROW_MISC_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <fast_matmult/common_asm.h>
#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

// ========================================================================

void serial_matmult(unsigned int M, unsigned int K, unsigned int N,
                    cblas_float *A, cblas_float *B, cblas_float *C);

void serial_matmult_sse(unsigned int M, unsigned int K, unsigned int N,
                        cblas_float *A, cblas_float *B, cblas_float *C);

// ========================================================================

void matmult_s_tiling_sse2(unsigned int M, unsigned int K, unsigned int N,
                           cblas_float *A, cblas_float *B, cblas_float *C);

void matmult_s_tiling_sse2_4x2(unsigned int M, unsigned int K, unsigned int N,
                               cblas_float *A, cblas_float *B, cblas_float *C);

// ========================================================================

void matrix_fast_matmult(unsigned int M, unsigned int K, unsigned int N,
                         cblas_float *A, cblas_float *B, cblas_float *C);

void matrix_fast_matmult_sse2_4x2(unsigned int M, unsigned int K, unsigned int N,
                                  cblas_float *A, cblas_float *B, cblas_float *C);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_ROW_MISC_H_ */
