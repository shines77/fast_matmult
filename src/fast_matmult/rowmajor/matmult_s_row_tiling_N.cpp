
#include <fast_matmult/rowmajor/matmult_s_row_tiling_N.h>

void matmult_s_row_tiling_MxK_N(const int M, const int N, const int K,
                                const cblas_float alpha,
                                const cblas_float *A, const int lda,
                                const cblas_float *B, const int ldb,
                                const cblas_float beta,
                                cblas_float *C, const int ldc)
{
    int m, n, k;
    int m_start, m_end;
    int n_start, n_end;
    int k_start, k_end;
    int m_max, n_max, k_max;
    int m_step, n_step, k_step;

    const cblas_float *A_, *B_;
    cblas_float *C_;
    cblas_float A_m_k;

    m_step = 4;
    k_step = 128;
    n_step = 512;

#if 0
    m_step = M;
    k_step = K;
    n_step = N;
#endif

    if (m_step > M)
        m_step = M;
    if (k_step > K)
        k_step = K;
    if (n_step > N)
        n_step = N;

    // 分块外层循环顺序: K, M, N
    TILING_OUTER_LOOP_BEGIN(k, m, n);
    //TILING_OUTER_LOOP_BEGIN(m, k, n);

        // 内层循环顺序: m, k - n
        TILING_INNER_LOOP_BEGIN(m, k);

        do {
            A_ = &A[(m + 1) * lda + k];
            A_m_k = A[m * lda + k];

            B_ = &B[k * ldb + n_start];
            C_ = &C[m * ldc + n_start];

            // 最内层循环: n
            for (n = n_start; n < n_end; ++n) {
                // C[m, n] += A[m, k] * B[k, n];
                //C[m * ldc + n] += A_m_k * B[k * ldb + n];
                (*C_++) += A_m_k * (*B_++);
            }
        } while (0);

        TILING_INNER_LOOP_END();

    TILING_OUTER_LOOP_END(k, m, n);
    //TILING_OUTER_LOOP_END(m, k, n);
}

void matmult_s_row_tiling_KxM_N(const int M, const int N, const int K,
                                const cblas_float alpha,
                                const cblas_float *A, const int lda,
                                const cblas_float *B, const int ldb,
                                const cblas_float beta,
                                cblas_float *C, const int ldc)
{
    int m, n, k;
    int m_start, m_end;
    int n_start, n_end;
    int k_start, k_end;
    int m_max, n_max, k_max;
    int m_step, n_step, k_step;

    const cblas_float *A_, *B_;
    cblas_float *C_;
    cblas_float A_m_k;

    m_step = 4;
    k_step = 128;
    n_step = 512;

#if 0
    m_step = M;
    k_step = K;
    n_step = N;
#endif

    if (m_step > M)
        m_step = M;
    if (k_step > K)
        k_step = K;
    if (n_step > N)
        n_step = N;

    // 分块外层循环顺序: K, M, N
    TILING_OUTER_LOOP_BEGIN(k, m, n);
    //TILING_OUTER_LOOP_BEGIN(m, k, n);

        // 内层循环顺序: k, m - n
        TILING_INNER_LOOP_BEGIN(k, m);

        do {
            A_ = &A[(m + 1) * lda + k];
            A_m_k = A[m * lda + k];

            B_ = &B[k * ldb + n_start];
            C_ = &C[m * ldc + n_start];

            // 最内层循环: n
            for (n = n_start; n < n_end; ++n) {
                // C[m, n] += A[m, k] * B[k, n];
                //C[m * ldc + n] += A_m_k * B[k * ldb + n];
                (*C_++) += A_m_k * (*B_++);
            }
        } while (0);

        TILING_INNER_LOOP_END();

    TILING_OUTER_LOOP_END(k, m, n);
    //TILING_OUTER_LOOP_END(m, k, n);
}

void matmult_s_row_tiling_k_n_m_KxM_N(const int M, const int N, const int K,
                                      const cblas_float alpha,
                                      const cblas_float *A, const int lda,
                                      const cblas_float *B, const int ldb,
                                      const cblas_float beta,
                                      cblas_float *C, const int ldc)
{
    int m, n, k;
    int m_start, m_end;
    int n_start, n_end;
    int k_start, k_end;
    int m_max, n_max, k_max;
    int m_step, n_step, k_step;

    const cblas_float *A_, *B_;
    cblas_float *C_;
    cblas_float A_m_k;

    m_step = 4;
    k_step = 128;
    n_step = 512;

#if 0
    m_step = M;
    k_step = K;
    n_step = N;
#endif

    if (m_step > M)
        m_step = M;
    if (k_step > K)
        k_step = K;
    if (n_step > N)
        n_step = N;

    // 分块外层循环顺序: K, N, M
    TILING_OUTER_LOOP_BEGIN(k, n, m);

        // 内层循环顺序: k, m - n
        TILING_INNER_LOOP_BEGIN(k, m);

        do {
            A_ = &A[(m + 1) * lda + k];
            A_m_k = A[m * lda + k];

            B_ = &B[k * ldb + n_start];
            C_ = &C[m * ldc + n_start];

            // 最内层循环: n
            for (n = n_start; n < n_end; ++n) {
                // C[m, n] += A[m, k] * B[k, n];
                //C[m * ldc + n] += A_m_k * B[k * ldb + n];
                (*C_++) += A_m_k * (*B_++);
            }
        } while (0);

        TILING_INNER_LOOP_END();

    TILING_OUTER_LOOP_END(k, n, m);
}
