
#include "fast_matmult/rowmajor/matmult_s_row_X.h"

/* 循环顺序: m, k - n */
void matmult_s_row_MxK_N(const int M, const int N, const int K,
                         const float_t alpha,
                         const float_t *A, const int lda,
                         const float_t *B, const int ldb,
                         const float_t beta,
                         float_t *C, const int ldc)
{
    float_t A_m_k;
    int m, n, k;

    for (m = 0; m < M; ++m) {
        for (k = 0; k < K; ++k) {
            // A_m_k = A[m, k];
            A_m_k = A[m * K + k];
            for (n = 0; n < N; ++n) {
                // C[m, n] += A_m_k * B[k, n];
                C[m * N + n] += A_m_k * B[k * N + n];
            }
        }
    }
}

/* 循环顺序: k, m - n */
void matmult_s_row_KxM_N(const int M, const int N, const int K,
                         const float_t alpha,
                         const float_t *A, const int lda,
                         const float_t *B, const int ldb,
                         const float_t beta,
                         float_t *C, const int ldc)
{
    float_t A_m_k;
    int m, n, k;

    for (k = 0; k < K; ++k) {
        for (m = 0; m < M; ++m) {
            // A_m_k = A[m, k];
            A_m_k = A[m * K + k];
            for (n = 0; n < N; ++n) {
                // C[m, n] += A_m_k * B[k, n];
                C[m * N + n] += A_m_k * B[k * N + n];
            }
        }
    }
}

/* 循环顺序: m, n - k */
void matmult_s_row_MxN_K(const int M, const int N, const int K,
                         const float_t alpha,
                         const float_t *A, const int lda,
                         const float_t *B, const int ldb,
                         const float_t beta,
                         float_t *C, const int ldc)
{
    float_t C_m_n;
    int m, n, k;

    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            C_m_n = (float_t)0.0;
            for (k = 0; k < K; ++k) {
                // C[m, n] += A[m, k] * B[k, n];
                C_m_n += A[m * K + k] * B[k * N + n];
            }
            C[m * N + n] = C_m_n;
        }
    }
}

/* 循环顺序: n, m - k */
void matmult_s_row_NxM_K(const int M, const int N, const int K,
                         const float_t alpha,
                         const float_t *A, const int lda,
                         const float_t *B, const int ldb,
                         const float_t beta,
                         float_t *C, const int ldc)
{
    float_t C_m_n;
    int m, n, k;

    for (n = 0; n < N; ++n) {
        for (m = 0; m < M; ++m) {
            C_m_n = (float_t)0.0;
            for (k = 0; k < K; ++k) {
                // C[m, n] += A[m, k] * B[k, n];
                C_m_n += A[m * K + k] * B[k * N + n];
            }
            C[m * N + n] = C_m_n;
        }
    }
}

/* 循环顺序: k, n - m */
void matmult_s_row_KxN_M(const int M, const int N, const int K,
                         const float_t alpha,
                         const float_t *A, const int lda,
                         const float_t *B, const int ldb,
                         const float_t beta,
                         float_t *C, const int ldc)
{
    float_t B_k_n;
    int m, n, k;

    for (k = 0; k < K; ++k) {
        for (n = 0; n < N; ++n) {
            // B_k_n = B[k, n];
            B_k_n = B[k * N + n];
            for (m = 0; m < M; ++m) {
                // C[m, n] += A[m, k] * B_k_n;
                C[m * N + n] += A[m * K + k] * B_k_n;
            }
        }
    }
}

/* 循环顺序: n, k - m */
void matmult_s_row_NxK_M(const int M, const int N, const int K,
                         const float_t alpha,
                         const float_t *A, const int lda,
                         const float_t *B, const int ldb,
                         const float_t beta,
                         float_t *C, const int ldc)
{
    float_t B_k_n;
    int m, n, k;

    for (n = 0; n < N; ++n) {
        for (k = 0; k < K; ++k) {
            // B_k_n = B[k, n];
            B_k_n = B[k * N + n];
            for (m = 0; m < M; ++m) {
                // C[m, n] += A[m, k] * B_k_n;
                C[m * N + n] += A[m * K + k] * B_k_n;
            }
        }
    }
}

/* 循环顺序: m, n - k */
void matmult_s_row_MxN_K_transB(const int M, const int N, const int K,
                                const float_t alpha,
                                const float_t *A, const int lda,
                                const float_t *B, const int ldb,
                                const float_t beta,
                                float_t *C, const int ldc)
{
    float_t C_m_n;
    int m, n, k;

    // 先转置矩阵B
    matrix_fast_transpose_NxN((float_t *)B, K, N);

    // 循环顺序: m, n - k
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            C_m_n = (float_t)0.0;
            for (k = 0; k < K; ++k) {
                // C[m, n] += A[m, k] * B'[k, n];
                C_m_n += A[m * K + k] * B[n * K + k];
            }
            C[m * N + n] = C_m_n;
        }
    }

    // 计算完以后再转置(还原)矩阵B
    matrix_fast_transpose_NxN((float_t *)B, N, K);
}

/* 循环顺序: n, m - k */
void matmult_s_row_NxM_K_transB(const int M, const int N, const int K,
                                const float_t alpha,
                                const float_t *A, const int lda,
                                const float_t *B, const int ldb,
                                const float_t beta,
                                float_t *C, const int ldc)
{
    float_t C_m_n;
    int m, n, k;

    // 先转置矩阵B
    matrix_fast_transpose_NxN((float_t *)B, K, N);

    // 循环顺序: n, m - k
    for (n = 0; n < N; ++n) {
        for (m = 0; m < M; ++m) {
            C_m_n = (float_t)0.0;
            for (k = 0; k < K; ++k) {
                // C[m, n] += A[m, k] * B'[k, n];
                C_m_n += A[m * K + k] * B[n * K + k];
            }
            C[m * N + n] = C_m_n;
        }
    }

    // 计算完以后再转置(还原)矩阵B
    matrix_fast_transpose_NxN((float_t *)B, N, K);
}
