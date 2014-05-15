
#include "fast_matmult/serial_matmult_colmajor_X.h"

namespace annlab {

/* ѭ��˳��: m, k - n */
void serial_matmult_colmajor_N_MxK(const int M, const int N, const int K,
                                   const float_t alpha,
                                   const float_t *A, const int lda,
                                   const float_t *B, const int ldb,
                                   const float_t beta,
                                   float_t *C, const int ldc)
{
    float_t A_m_k;
    int m, n, k;
    // m��kѭ��˳��ɽ���
    for (m = 0; m < M; ++m) {
        for (k = 0; k < K; ++k) {
            // A_m_k = A[m, k];
            A_m_k = A[k * M + m];
            for (n = 0; n < N; ++n) {
                // C[m, n] += A_m_k * B[k, n];
                C[n * M + m] += A_m_k * B[n * K + k];
            }
        }
    }
}

/* ѭ��˳��: k, m - n */
void serial_matmult_colmajor_N_KxM(const int M, const int N, const int K,
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
            A_m_k = A[k * M + m];
            for (n = 0; n < N; ++n) {
                // C[m, n] += A_m_k * B[k, n];
                C[n * M + m] += A_m_k * B[n * K + k];
            }
        }
    }
}

/* ѭ��˳��: m, n - k */
void serial_matmult_colmajor_K_MxN(const int M, const int N, const int K,
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
                C_m_n += A[k * M + m] * B[n * K + k];
            }
            C[n * M + m] = C_m_n;
        }
    }
}

/* ѭ��˳��: n, m - k */
void serial_matmult_colmajor_K_NxM(const int M, const int N, const int K,
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
                C_m_n += A[k * M + m] * B[n * K + k];
            }
            C[n * M + m] = C_m_n;
        }
    }
}

/* ѭ��˳��: k, n - m */
void serial_matmult_colmajor_M_KxN(const int M, const int N, const int K,
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
            B_k_n = B[n * K + k];
            for (m = 0; m < M; ++m) {
                // C[m, n] += A[m, k] * B_k_n;
                C[n * M + m] += A[k * M + m] * B_k_n;
            }
        }
    }
}

/* ѭ��˳��: n, k - m */
void serial_matmult_colmajor_M_NxK(const int M, const int N, const int K,
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
            B_k_n = B[n * K + k];
            for (m = 0; m < M; ++m) {
                // C[m, n] += A[m, k] * B_k_n;
                C[n * M + m] += A[k * M + m] * B_k_n;
            }
        }
    }
}

/* ѭ��˳��: m, n - k */
void serial_matmult_colmajor_K_MxN_transA(const int M, const int N, const int K,
                                          const float_t alpha,
                                          const float_t *A, const int lda,
                                          const float_t *B, const int ldb,
                                          const float_t beta,
                                          float_t *C, const int ldc)
{
    float_t C_m_n;
    int m, n, k;

    // ��ת�þ���A
    matrix_fast_transpose_NxN((float_t *)A, M, K);

    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            C_m_n = (float_t)0.0;
            for (k = 0; k < K; ++k) {
                // C[m, n] += A'[m, k] * B[k, n];
                C_m_n += A[m * K + k] * B[n * K + k];
            }
            C[n * M + m] = C_m_n;
        }
    }

    // �������Ժ���ת��(��ԭ)����A
    matrix_fast_transpose_NxN((float_t *)A, K, M);
}

/* ѭ��˳��: n, m - k */
void serial_matmult_colmajor_K_NxM_transA(const int M, const int N, const int K,
                                          const float_t alpha,
                                          const float_t *A, const int lda,
                                          const float_t *B, const int ldb,
                                          const float_t beta,
                                          float_t *C, const int ldc)
{
    float_t C_m_n;
    int m, n, k;

    // ��ת�þ���A
    matrix_fast_transpose_NxN((float_t *)A, M, K);

    for (n = 0; n < N; ++n) {
        for (m = 0; m < M; ++m) {
            C_m_n = (float_t)0.0;
            for (k = 0; k < K; ++k) {
                // C[m, n] += A'[m, k] * B[k, n];
                C_m_n += A[m * K + k] * B[n * K + k];
            }
            C[n * M + m] = C_m_n;
        }
    }

    // �������Ժ���ת��(��ԭ)����A
    matrix_fast_transpose_NxN((float_t *)A, K, M);
}

}  /* namespace annlab */
