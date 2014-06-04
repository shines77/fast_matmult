
#include <fast_matmult/rowmajor/matmult_s_row_tiling_K.h>

void matmult_s_row_tiling_NxM_K(const int M, const int N, const int K,
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
    //float_t *C_;
    cblas_float C_m_n;

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
    //TILING_OUTER_LOOP_BEGIN(n, m, k);

        // 内层循环顺序: n, m - k
        TILING_INNER_LOOP_BEGIN(n, m);

        do {
            A_ = &A[m * lda + k_start];

            B_ = &B[k_start * ldb + n];
            //C_ = &C[m * ldc + n];

            C_m_n = (cblas_float)0.0;

            // 最内层循环: k
            for (k = k_start; k < k_end; ++k) {
                // C[m, n] += A[m, k] * B[k, n];
                //C_m_n += A[m * K + k] * B[k * N + n];
                C_m_n += (*A_++) * (*B_);
                B_ += ldb;
            }
            C[m * ldc + n] += C_m_n;
        } while (0);

        TILING_INNER_LOOP_END();

    TILING_OUTER_LOOP_END(k, m, n);
    //TILING_OUTER_LOOP_END(n, m, k);
}

void matmult_s_row_tiling_MxN_K_transB(unsigned int M, unsigned int K, unsigned int N,
                                       cblas_float *A, cblas_float *B, cblas_float *C)
{
    unsigned int m, n, k;
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int k_start, k_end;
    unsigned int m_step, n_step, k_step;
    cblas_float *C_ = NULL, *B_, *A_;
    cblas_float C_m_n;

    //
    // matrix multiplication: C1 = A * B
    //
    // [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
    // [D][E][F] * [I][J]  =  [D*G + E*I + F*K][D*H + E*J + F*L]
    //             [K][L]
    //

    // 先转置矩阵B
    matrix_fast_transpose_NxN(B, K, N);

    m_step = 4;
    k_step = 512;
    n_step = 128;

#if 0
    m_step = 2;
    k_step = 1024;
    n_step = 64;
#endif

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

#if 1
    // 分块循环顺序N, M, K, 内层循环顺序m, n, k
    n_start = 0;
    n_end   = n_start + n_step;
    if (n_end > N)
        n_end = N;
    while (n_start < N) {
        m_start = 0;
        m_end   = m_start + m_step;
        if (m_end > M)
            m_end = M;
        while (m_start < M) {
            k_start = 0;
            k_end   = k_start + k_step;
            if (k_end > K)
                k_end = K;
            while (k_start < K) {
                // m_step (0 - 4),
                // C[m, n] need m_step * ((n_step * 8) / 4096) = 4 x 1 TLB entries,
                // A[m, k] only need m_step x |k_step / nPageSize| = 4 x 1 TLB entries.
                for (m = m_start; m < m_end; ++m) {
                    // n_step (0 - 128),
                    // (nPageSize = 4096 bytes)
                    // B[n, k] need n_step * ((k_step * 8) / nPageSize) = n_step x 1 = 128 x 1 TLB entries,
                    // C[m, n] need m_step * ((n_step * 8) / nPageSize) = m_step x 1 = 4 x 1 TLB entries.
                    for (n = n_start; n < n_end; ++n) {
    #if 1
                        A_ = &A[m * K + k_start];
                        B_ = &B[n * K + k_start];
                        //C_ = &C[m * N + n];
                        /*
                        if (k_step >= K)
                            C_m_n = (float_t)0.0;
                        else
                        //*/
                            C_m_n = C[m * N + n];
                        //_mm_prefetch((char *)A_, _MM_HINT_NTA);
                        //_mm_prefetch((char *)B_, _MM_HINT_NTA);

                        // k_step (0 - 512),
                        // B[n, k] need 128 x 1 TLB entries,
                        // A[m, k] need m_step x 1 TLB entries if k_step < nPageSize.
                        for (k = k_start; k < k_end; ++k) {
                            // C[m, n] += A[m, k] * B'[n, k];
                            //C[m * N + n] += A[m * K + k] * B[n * K + k];
                            C_m_n += (*A_++) * (*B_++);
                        }
                        C[m * N + n] = C_m_n;
    #elif 0
                        // SSE2 汇编
    #endif
                    }
                }
                k_start += k_step;
                k_end   += k_step;
                if (k_end > K)
                    k_end = K;
            }
            m_start += m_step;
            m_end   += m_step;
            if (m_end > M)
                m_end = M;
        }
        n_start += n_step;
        n_end   += n_step;
        if (n_end > N)
            n_end = N;
    }
#elif 0
    // 分块循环顺序M, N, K, 内层循环顺序m, n, k
    m_start = 0;
    m_end   = m_start + m_step;
    if (m_end > M)
        m_end = M;
    while (m_start < M) {
        n_start = 0;
        n_end   = n_start + n_step;
        if (n_end > N)
            n_end = N;
        while (n_start < N) {
            k_start = 0;
            k_end   = k_start + k_step;
            if (k_end > K)
                k_end = K;
            while (k_start < K) {
                // m_step (0 - 4),
                // C[m, n] need m_step * ((n_step * 8) / 4096) = 4 x 1 TLB entries,
                // A[m, k] only need m_step x |k_step / nPageSize| = 4 x 1 TLB entries.
                for (m = m_start; m < m_end; ++m) {
                //for (n = n_start; n < n_end; ++n) {
                    // n_step (0 - 128),
                    // (nPageSize = 4096 bytes)
                    // B[n, k] need n_step * ((k_step * 8) / nPageSize) = n_step x 1 = 128 x 1 TLB entries,
                    // C[m, n] need m_step * ((n_step * 8) / nPageSize) = m_step x 1 = 4 x 1 TLB entries.
                    for (n = n_start; n < n_end; ++n) {
                    //for (m = m_start; m < m_end; ++m) {
    #if 1
                        A_ = &A[m * K + k_start];
                        B_ = &B[n * K + k_start];
                        //C_ = &C[m * N + n];
                        //C_m_n = (float_t)0.0;
                        C_m_n = C[m * N + n];
                        //_mm_prefetch((char *)A_, _MM_HINT_NTA);
                        //_mm_prefetch((char *)B_, _MM_HINT_NTA);

                        // k_step (0 - 512),
                        // B[n, k] need 128 x 1 TLB entries,
                        // A[m, k] need m_step x 1 TLB entries if k_step < nPageSize.
                        for (k = k_start; k < k_end; ++k) {
                            // C[m, n] += A[m, k] * B'[n, k];
                            //C[m * N + n] += A[m * K + k] * B[n * K + k];
                            C_m_n += (*A_++) * (*B_++);
                        }
                        C[m * N + n] = C_m_n;
    #elif 0
                        // SSE2 汇编
    #endif
                    }
                }
                k_start += k_step;
                k_end   += k_step;
                if (k_end > K)
                    k_end = K;
            }
            n_start += n_step;
            n_end   += n_step;
            if (n_end > N)
                n_end = N;
        }
        m_start += m_step;
        m_end   += m_step;
        if (m_end > M)
            m_end = M;
    }
#else
    // 分块循环顺序M, N, K, 内层循环顺序m, n, k
    k_start = 0;
    k_end   = k_start + k_step;
    if (k_end > K)
        k_end = K;
    while (k_start < K) {
        m_start = 0;
        m_end   = m_start + m_step;
        if (m_end > M)
            m_end = M;
        while (m_start < M) {
            n_start = 0;
            n_end   = n_start + n_step;
            if (n_end > N)
                n_end = N;
            while (n_start < N) {
                // m_step (0 - 4),
                // C[m, n] need m_step * ((n_step * 8) / 4096) = 4 x 1 TLB entries,
                // A[m, k] only need m_step x |k_step / nPageSize| = 4 x 1 TLB entries.
                for (m = m_start; m < m_end; ++m) {
                    // n_step (0 - 128),
                    // (nPageSize = 4096 bytes)
                    // B[n, k] need n_step * ((k_step * 8) / nPageSize) = n_step x 1 = 128 x 1 TLB entries,
                    // C[m, n] need m_step * ((n_step * 8) / nPageSize) = m_step x 1 = 4 x 1 TLB entries.
                    for (n = n_start; n < n_end; ++n) {
    #if 1
                        A_ = &A[m * K + k_start];
                        B_ = &B[n * K + k_start];
                        //C_ = &C[m * N + n];
                        //C_m_n = (float_t)0.0;
                        C_m_n = C[m * N + n];
                        //_mm_prefetch((char *)A_, _MM_HINT_NTA);
                        //_mm_prefetch((char *)B_, _MM_HINT_NTA);

                        // k_step (0 - 512),
                        // B[n, k] need 128 x 1 TLB entries,
                        // A[m, k] need m_step x 1 TLB entries if k_step < nPageSize.
                        for (k = k_start; k < k_end; ++k) {
                            // C[m, n] += A[m, k] * B'[n, k];
                            //C[m * N + n] += A[m * K + k] * B[n * K + k];
                            C_m_n += (*A_++) * (*B_++);
                        }
                        C[m * N + n] = C_m_n;
    #elif 0
                        // SSE2 汇编
    #endif
                    }
                }
                n_start += n_step;
                n_end   += n_step;
                if (n_end > N)
                    n_end = N;
            }
            m_start += m_step;
            m_end   += m_step;
            if (m_end > M)
                m_end = M;
        }
        k_start += k_step;
        k_end   += k_step;
        if (k_end > K)
            k_end = K;
    }
#endif

    // 计算完以后再转置(还原)矩阵B
    matrix_fast_transpose_NxN(B, N, K);
}
