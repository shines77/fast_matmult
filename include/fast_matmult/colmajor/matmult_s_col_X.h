
#ifndef _MATMULT_S_COL_X_H_
#define _MATMULT_S_COL_X_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#if (_WIN32 | _WIN64) && _MSC_VER
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

#include <fast_matmult/fast_matmult.h>

#ifdef __cplusplus
extern "C" {
#endif

// ========================================================================

void matmult_s_col_MxN_K(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

void matmult_s_col_NxM_K(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

void matmult_s_col_MxK_N(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

void matmult_s_col_KxM_N(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

void matmult_s_col_NxK_M(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

void matmult_s_col_KxN_M(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

// ========================================================================

void matmult_s_col_MxN_K_transA(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

void matmult_s_col_NxM_K_transA(
    const int M, const int N, const int K,
    const cblas_float alpha,
    const cblas_float *A, const int lda,
    const cblas_float *B, const int ldb,
    const cblas_float beta,
    cblas_float *C, const int ldc);

// ========================================================================

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_S_COL_X_H_ */
