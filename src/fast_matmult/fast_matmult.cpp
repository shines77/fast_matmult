
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _MSC_VER
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#else
// non-msvc xmm & emm define header files
// TODO: xmm and emm header files
#endif

#include "fast_matmult/fast_matmult.h"
#include "fast_matmult/common_asm.h"
#include "fast_matmult/stop_watch.h"
#include "fast_matmult/huge_tlb.h"

#include "fast_matmult/matmult_s_row_X.h"
#include "fast_matmult/matmult_s_col_X.h"
#include "fast_matmult/matmult_s_row_tiling_K.h"
#include "fast_matmult/matmult_s_row_tiling_N.h"

using namespace annlab;

#define DISPLAY_MATRIX_COMPARE      0

#define DOUBLE_NUMS_PER_LOOP        16
#define MM_PREFETCH_OFFSET          1024
#define MM_PREFETCH_OFFSET_V        (MM_PREFETCH_OFFSET / sizeof(double))
#define USE_SEE2_WRITE_PREFRTCH     0

#define MM_PREFETCH_OFFSET2         128
#define MM_PREFETCH_OFFSET_V2       (MM_PREFETCH_OFFSET2 / sizeof(double))

/* 0=ASM嵌入SSE2代码, 1=纯C代码, 2=编译器级别SSE2代码, 3=其他(默认) */
#if (_WIN32 | _WIN64) && _MSC_VER
#  define _MULT_SSE2_MODE           2
#else
#  define _MULT_SSE2_MODE           1
#endif

#define MIN(a, b)       (((a) < (b)) ? (a) : (b))
#define MAX(a, b)       (((a) > (b)) ? (a) : (b))

float_t *matrix_malloc(unsigned int M, unsigned int N,
                       unsigned int alignment)
{
    return (float_t *)::_aligned_malloc(M * N * sizeof(float_t), alignment);
}

float_t *matrix_malloc_ex(unsigned int M, unsigned int N,
                          unsigned int alignment,
                          eMatInitFcn initFcn /* = MatInitZeros */,
                          float_t fillValue /* = 0.0 */,
                          eMatrixItemOrder order /* = MatItemOrderAsc */)
{
    float_t *A;
    A = (float_t *)::_aligned_malloc(M * N * sizeof(float_t), alignment);
    if (A != NULL) {
        // init elements
        matrix_init_elements(A, M, N, initFcn, fillValue, order);
    }
    return A;
}

void matrix_free(float_t *A)
{
    ::_aligned_free((void *)A);
}

void matrix_init_elements(float_t *A, unsigned int M, unsigned int N,
                          eMatInitFcn initFcn /* = MatInitZeros */,
                          float_t fillValue /* = 0.0 */,
                          eMatrixItemOrder order /* = MatItemOrderAsc */)
{
    unsigned int m, n;
    int index, value;

    if (initFcn == MatInitZeros) {
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                index = n * M + m;
                A[index] = (float_t)0.0;
            }
        }
    }
    else if (initFcn == MatInitOnes) {
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                index = n * M + m;
                A[index] = (float_t)1.0;
            }
        }
    }
    else if (initFcn == MatInitRands) {
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                index = n * M + m;
                A[index] = (float_t)(::rand() / (float_t)RAND_MAX * (float_t)2.0 - (float_t)1.0);
            }
        }
    }
    else if (initFcn == MatInitRands_Positive) {
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                index = n * M + m;
                A[index] = (float_t)(::rand() / (float_t)RAND_MAX);
            }
        }
    }
    else if (initFcn == MatInitOrder) {
        if (order == MatItemOrderAsc) {
            for (n = 0; n < N; ++n) {
                for (m = 0; m < M; ++m) {
                    index = n * M + m;
                    A[index] = (float_t)index;
                }
            }
        }
        else if (order == MatItemOrderDesc) {
            for (n = 0; n < N; ++n) {
                for (m = 0; m < M; ++m) {
                    index = n * M + m;
                    value = m * N + n;
                    A[index] = (float_t)value;
                }
            }
        }
        else {
            // unknown order type
        }
    }
    else if (initFcn == MatInitSpecify) {
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                index = n * M + m;
                A[index] = (float_t)fillValue;
            }
        }
    }
    else {
        // unknown init function
    }
}

bool matrix_compare(const float_t *A, const float_t *B, unsigned int M, unsigned int N,
                    int *diff_nums /* = NULL */,
                    eMatrixItemOrder order /* = MatItemOrderAsc */)
{
    bool cmp_result = false;
    const double err_epsilon = FLOAT_T_EPSINON_TEST;
    unsigned int m, n;
    float_t cell_val1, cell_val2;
    int diff_cnt;
    if (order == MatItemOrderAsc) {
        diff_cnt = 0;
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                cell_val1 = A[n * M + m];
                cell_val2 = B[n * M + m];
                if (::fabs(cell_val1 - cell_val2) > err_epsilon)
                    diff_cnt++;
            }
        }
        if (diff_cnt == 0)
            cmp_result = true;
    }
    else if (order == MatItemOrderTransAsc) {
        diff_cnt = 0;
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                cell_val1 = A[n * M + m];
                cell_val2 = B[m * N + n];
                if (::fabs(cell_val1 - cell_val2) > err_epsilon)
                    diff_cnt++;
            }
        }
        if (diff_cnt == 0)
            cmp_result = true;
    }

    if (diff_nums != NULL)
        *diff_nums = diff_cnt;
    return cmp_result;
}

bool matrix_transpose_verify(float_t *A, unsigned int M, unsigned int N,
                             int *err_nums /* = NULL */,
                             eMatrixItemOrder order /* = MatItemOrderAsc */)
{
    bool verify_ok = false;
    const double err_epsilon = FLOAT_T_EPSILON;
    unsigned int m, n;
    float_t cell_val;
    int index, err_cnt;

    if (order == MatItemOrderAsc) {
        err_cnt = 0;
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                cell_val = A[n * M + m];
                index = n * M + m;
                if (::fabs(cell_val - (float_t)index) > err_epsilon)
                    err_cnt++;
            }
        }
        if (err_cnt == 0)
            verify_ok = true;
    }
    else if (order == MatItemOrderDesc) {
        err_cnt = 0;
        for (n = 0; n < N; ++n) {
            for (m = 0; m < M; ++m) {
                cell_val = A[n * M + m];
                index = m * N + n;
                if (::fabs(cell_val - (float_t)index) > err_epsilon)
                    err_cnt++;
            }
        }
        if (err_cnt == 0)
            verify_ok = true;
    }
    else {
        // unknown order type
    }

    if (err_nums != NULL)
        *err_nums = err_cnt;
    return verify_ok;
}

void matrix_fast_transpose_NxN(float_t *A, unsigned int M, unsigned int N)
{
    unsigned int m, n;
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int m_step, n_step;
    //unsigned int n_start_min = 0, n_start_max = 0;
    double temp1, temp2;

    m_step = 128;
    n_step = 128;

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
            // (0 - 128), 128 x 8 = 1024 bytes(page size), need 1 TLB entry
            for (m = m_start; m < MIN(m_end, M); ++m) {
                // (0 - 128), need 128 TLB entry
                if (MAX(n_start, m + 1) > MIN(n_end, N))
                //if (n_start_min > n_start_max)
                    break;
                //n_start_min = MAX(n_start, m + 1);
                //n_start_max = MIN(n_end, N);
                //for (n = n_start_min; n < n_start_max; ++n) {
                for (n = MAX(n_start, m + 1); n < MIN(n_end, N); ++n) {
#if 0
                    temp1 = A[n * M + m];
                    A[n * M + m] = A[m * N + n];
                    A[m * N + n] = temp1;
#else
                    temp1 = A[m * N + n];
                    temp2 = A[n * M + m];
                    A[n * M + m] = temp1;
                    A[m * N + n] = temp2;
#endif
                }
            }
            n_start += n_step;
            n_end   += n_step;
        }
        m_start += m_step;
        m_end   += m_step;
    }
}

#define MATRIX_AB_NEED_TRANSPOSE     0

void matrix_fast_matmult(unsigned int M, unsigned int K, unsigned int N,
                         float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int k_start, k_end;
    unsigned int m_step, n_step, k_step;
    double temp = 0.0, temp2 = 0.0;
    double elapsedTime1 = 0.0;
    int err_nums = 0;
    bool verify_ok = false;

    //
    // Step1: if need, first transpose the matrix B
    //

#if MATRIX_AB_NEED_TRANSPOSE

    stop_watch stopWatch;

    m_step = 128;
    n_step = 128;

    stopWatch.start();

#if 0
    m_start = 0;
    m_end   = m_start + m_step;

    while (m_start <= M) {
        n_start = 0;
        n_end   = n_start + n_step;
        while (n_start <= N) {
            // (0 - 128), 128 x 8 = 1024 bytes(1/4 page size), only need one TLB entry
            for (m = m_start; m < MIN(m_end, M); ++m) {
                // (0 - 128), need 128 TLB entry
                if (MAX(n_start, m + 1) > MIN(n_end, N))
                //if (n_start_min > n_start_max)
                    break;
                //n_start_min = MAX(n_start, m + 1);
                //n_start_max = MIN(n_end, N);
                //for (n = n_start_min; n < n_start_max; ++n) {
                for (n = MAX(n_start, m + 1); n < MIN(n_end, N); ++n) {
#if 0
                    temp = B[n * M + m];
                    B[n * M + m] = B[m * N + n];
                    B[m * N + n] = temp;
#else
                    temp  = B[m * N + n];
                    temp2 = B[n * M + m];
                    B[n * M + m] = temp;
                    B[m * N + n] = temp2;
#endif
                }
            }
            n_start += n_step;
            n_end   += n_step;
        }
        m_start += m_step;
        m_end   += m_step;
    }
#else
    n_start = 0;
    n_end   = n_start + n_step;

    while (n_start <= N) {
        m_start = 0;
        m_end   = m_start + m_step;
        while (m_start <= M) {
            // (0 - 128), need 128 TLB entry
            for (n = n_start; n < MIN(n_end, N); ++n) {
                if (MAX(m_start, n + 1) > MIN(m_end, M))
                    break;
                // (0 - 128), 128 x 8 = 1024 bytes(1/4 page size), only need one TLB entry
                for (m = MAX(m_start, n + 1); m < MIN(m_end, M); ++m) {
#if 0
                    temp = B[n * M + m];
                    B[n * M + m] = B[m * N + n];
                    B[m * N + n] = temp;
#else
                    temp  = B[n * M + m];
                    temp2 = B[m * N + n];
                    B[m * N + n] = temp;
                    B[n * M + m] = temp2;
#endif
                }
            }
            m_start += m_step;
            m_end   += m_step;
        }
        n_start += n_step;
        n_end   += n_step;
    }
#endif

    stopWatch.stop();
    elapsedTime1 = stopWatch.getMillisec();

    err_nums = 0;
    verify_ok = matrix_transpose_verify(B, K, N, &err_nums, MatItemOrderDesc);

#else  /* !MATRIX_AB_NEED_TRANSPOSE */

    elapsedTime1 = 0.0;
    verify_ok = true;

#endif  /* MATRIX_AB_NEED_TRANSPOSE */

    //matrix_init_elements(C, M, N, MatInitZeros);

    //
    // Step2: matrix multiplication, C1 = A * B
    //
    // [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
    // [D][E][F] * [I][J]  =  [D*G + E*I + F*K][D*H + E*J + F*L]
    //             [K][L]
    //

#if CBLAS_ROW_MAJOR
    float_t A_m_k;

    // (1.185 seconds)
    m_step = 32;
    k_step = 112;
    n_step = 1024;

    // (1.184-1.190 seconds)
    m_step = 64;
    k_step = 96;
    n_step = 1024;

    // (1.185-1.190 seconds)
    m_step = 32;
    k_step = 116;
    n_step = 1024;

    // (1.187-1.192 seconds)
    m_step = 128;
    k_step = 64;
    n_step = 1024;

    // (1.188-1.194 seconds)
    m_step = 64;
    k_step = 64;
    n_step = 1024;

    // (1.184-1.19 seconds)
    m_step = 64;
    k_step = 96;
    n_step = 1024;

    // (1.185-1.190 seconds)
    m_step = 32;
    k_step = 96;
    n_step = 1024;

    // (1.184-1.190 seconds)
    m_step = 1;
    k_step = 128;
    n_step = 1024;

    // (1.188-1.190 seconds)
    m_step = 16;
    k_step = 112;
    n_step = 1024;

    // (1.186-1.190 seconds)
    m_step = 1;
    k_step = 112;
    n_step = 1024;

    // (1.184-1.190 seconds)
    m_step = 1;
    k_step = 120;
    n_step = 1024;

    // (1.183-1.190 seconds)
    m_step = 2;
    k_step = 120;
    n_step = 1024;

    // (1.185-1.190 seconds), 分块循环顺序M, K, N
    m_step = 96;
    k_step = 32;
    n_step = 1024;

    m_step = 2;
    k_step = 120;
    n_step = 1024;

    m_step = 32;
    k_step = 96;
    n_step = 1024;

    // (1.080-1.090 seconds)
    m_step = 2;
    k_step = 256;
    n_step = 1024;

#if 0
    // (1.056 seconds)
    m_step = 4;
    k_step = 192;
    n_step = 512;
#endif

    // (1.065-1.072 seconds)
    m_step = 4;
    k_step = 224;
    n_step = 512;

    // (1.055-1.060 seconds)
    m_step = 4;
    k_step = 96;
    n_step = 512;

    // fastest (1.049-1.052 seconds)
    m_step = 4;
    k_step = 128;
    n_step = 512;

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

#if 1
    float_t *C_, *B_, *A_;
    // 分块循环顺序K, M, N
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
                // (0 - 32), need 32 TLB entry
                for (k = k_start; k < k_end; ++k) {
                    // (0 - 112), need n_step * 2(112 x 2) TLB entry
                    for (m = m_start; m < m_end; ++m) {
                        //
                        A_ = &A[(m + 1) * K + k];
                        _mm_prefetch((char *)A_, _MM_HINT_NTA);
#if MATRIX_AB_NEED_TRANSPOSE
                        A_m_k = A[k * M + m];
#else
                        A_m_k = A[m * K + k];
#endif  /* MATRIX_AB_NEED_TRANSPOSE */
                        if (fabs(A_m_k) > (float_t)FLOAT_T_EPSILON) {
                            C_ = &C[m * N + n_start];
                            B_ = &B[k * N + n_start];
#if (_MULT_SSE2_MODE == 1)
                            // (0 - 512), 512 x 8 = 4096 bytes(1 page size), only need 1 TLB entry,
                            // but in C[m, n], need n_step TLB entry too.
                            for (n = n_start; n < n_end; ++n) {
                                // C[m, n] += A[m, k] * B[k, n];
                                //C[m * N + n] += A_m_k * B[k * N + n];
                                (*C_++) += A_m_k * (*B_++);
                            }
#elif (_MULT_SSE2_MODE == 2)
                            __m128d tmp;
                            __m128d mm1, mm2, mm3, mm4;
    #if USE_SEE2_WRITE_PREFRTCH
                            __m128d mma, mmb;
                            //__m128d mmc, mmd;
    #endif
    #if (DOUBLE_NUMS_PER_LOOP == 16)
                            __m128d mm5, mm6, mm7, mm8;
    #endif
                            tmp = _mm_setzero_pd();
                            tmp = _mm_loadl_pd(tmp, &A_m_k);
                            tmp = _mm_shuffle_pd(tmp, tmp, 0);

                            for (n = (n_end - n_start) / DOUBLE_NUMS_PER_LOOP;
                                 n > 0;
                                 C_ += DOUBLE_NUMS_PER_LOOP, B_ += DOUBLE_NUMS_PER_LOOP, --n) {
    #if USE_SEE2_WRITE_PREFRTCH
                                // prefetchnta [addr + offset]
                                _mm_prefetch(((char *)B_ + MM_PREFETCH_OFFSET2), _MM_HINT_NTA);
                                _mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2), _MM_HINT_NTA);

                                mm1 = _mm_load_pd(B_   );
                                mm2 = _mm_load_pd(B_ + 2);

                                mma = _mm_load_pd(C_);
                                mmb = _mm_load_pd(C_ + 2);

                                mm1 = _mm_mul_pd(mm1, tmp);
                                mm2 = _mm_mul_pd(mm2, tmp);

                                mm3 = _mm_load_pd(B_ + 4);
                                mm4 = _mm_load_pd(B_ + 6);

                                mm3 = _mm_mul_pd(mm3, tmp);
                                mm4 = _mm_mul_pd(mm4, tmp);

                                //_mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2), _MM_HINT_NTA);

                                //_mm_prefetch(((char *)B_ + MM_PREFETCH_OFFSET2 + 32), _MM_HINT_NTA);
                                //_mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2 + 32), _MM_HINT_NTA);

                                mm1 = _mm_add_pd(mm1, mma);
                                mm2 = _mm_add_pd(mm2, mmb);

                                mm3 = _mm_add_pd(mm3, *((__m128d *)(C_ + 4)));
                                mm4 = _mm_add_pd(mm4, *((__m128d *)(C_ + 6)));

                                _mm_store_pd(C_,     mm1);
                                _mm_store_pd(C_ + 2, mm2);

        #if (DOUBLE_NUMS_PER_LOOP == 16)
                                _mm_prefetch(((char *)B_ + MM_PREFETCH_OFFSET2 + 64), _MM_HINT_NTA);
                                _mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2 + 64), _MM_HINT_NTA);

                                mm5 = _mm_load_pd(B_ + 8);
                                mm6 = _mm_load_pd(B_ + 10);
        #endif
                                _mm_store_pd(C_ + 4, mm3);
                                _mm_store_pd(C_ + 6, mm4);

        #if (DOUBLE_NUMS_PER_LOOP == 16)
                                mma = _mm_load_pd(C_ + 8);
                                mmb = _mm_load_pd(C_ + 10);

                                mm5 = _mm_mul_pd(mm5, tmp);
                                mm6 = _mm_mul_pd(mm6, tmp);

                                mm7 = _mm_load_pd(B_ + 12);
                                mm8 = _mm_load_pd(B_ + 14);

                                mm5 = _mm_add_pd(mm5, mma);
                                mm6 = _mm_add_pd(mm6, mmb);

                                //_mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2 + 64), _MM_HINT_NTA);

                                //_mm_prefetch(((char *)B_ + MM_PREFETCH_OFFSET2 + 96), _MM_HINT_NTA);
                                //_mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2 + 96), _MM_HINT_NTA);

                                mm7 = _mm_mul_pd(mm7, tmp);
                                mm8 = _mm_mul_pd(mm8, tmp);

                                _mm_store_pd(C_ + 8,  mm5);
                                _mm_store_pd(C_ + 10, mm6);

                                mm7 = _mm_add_pd(mm7, *((__m128d *)(C_ + 12)));
                                mm8 = _mm_add_pd(mm8, *((__m128d *)(C_ + 14)));

                                _mm_store_pd(C_ + 12, mm7);
                                _mm_store_pd(C_ + 14, mm8);
        #endif
    #else /* !USE_SEE2_WRITE_PREFRTCH */
                                // prefetchnta [addr + offset]
                                //_mm_prefetch(((char *)B_ + MM_PREFETCH_OFFSET2), _MM_HINT_NTA);
                                //_mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2), _MM_HINT_NTA);

                                mm1 = _mm_load_pd(B_   );
                                mm2 = _mm_load_pd(B_ + 2);

                                mm3 = _mm_load_pd(B_ + 4);
                                mm4 = _mm_load_pd(B_ + 6);

                                mm1 = _mm_mul_pd(mm1, tmp);
                                mm2 = _mm_mul_pd(mm2, tmp);

                                mm3 = _mm_mul_pd(mm3, tmp);
                                mm4 = _mm_mul_pd(mm4, tmp);

                                mm1 = _mm_add_pd(mm1, *((__m128d *)(C_))   );
                                mm2 = _mm_add_pd(mm2, *((__m128d *)(C_ + 2)));

                                mm3 = _mm_add_pd(mm3, *((__m128d *)(C_ + 4)));
                                mm4 = _mm_add_pd(mm4, *((__m128d *)(C_ + 6)));

                                _mm_store_pd(C_,     mm1);
                                _mm_store_pd(C_ + 2, mm2);

        #if (DOUBLE_NUMS_PER_LOOP == 16)
                                //_mm_prefetch(((char *)B_ + MM_PREFETCH_OFFSET2 + 64), _MM_HINT_NTA);
                                //_mm_prefetch(((char *)C_ + MM_PREFETCH_OFFSET2 + 64), _MM_HINT_NTA);

                                mm5 = _mm_load_pd(B_ + 8);
                                mm6 = _mm_load_pd(B_ + 10);
        #endif
                                _mm_store_pd(C_ + 4, mm3);
                                _mm_store_pd(C_ + 6, mm4);

        #if (DOUBLE_NUMS_PER_LOOP == 16)
                                mm5 = _mm_mul_pd(mm5, tmp);
                                mm6 = _mm_mul_pd(mm6, tmp);

                                mm7 = _mm_load_pd(B_ + 12);
                                mm8 = _mm_load_pd(B_ + 14);

                                mm5 = _mm_add_pd(mm5, *((__m128d *)(C_ +  8)));
                                mm6 = _mm_add_pd(mm6, *((__m128d *)(C_ + 10)));

                                mm7 = _mm_mul_pd(mm7, tmp);
                                mm8 = _mm_mul_pd(mm8, tmp);

                                _mm_store_pd(C_ + 8,  mm5);
                                _mm_store_pd(C_ + 10, mm6);

                                mm7 = _mm_add_pd(mm7, *((__m128d *)(C_ + 12)));
                                mm8 = _mm_add_pd(mm8, *((__m128d *)(C_ + 14)));

                                _mm_store_pd(C_ + 12, mm7);
                                _mm_store_pd(C_ + 14, mm8);
        #endif
                                /*
                                __asm {
                                    mov     eax, 64
L001:
                                    dec     eax
                                    //_emit   0x3E
                                    jne     L001
                                }
                                //*/
    #endif  /* USE_SEE2_WRITE_PREFRTCH */
                            }
#endif  /* _MULT_SSE2_MODE == 2 */
                        }
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
#elif 1
    float_t *C_, *B_;
    // 分块循环顺序K, M, N
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
            //n_end   = n_start + N;
            n_end   = n_start + n_step;
            if (n_end > N)
                n_end = N;
            while (n_start < N) {
                // (0 - 112), need n_step * 2(112 x 2) TLB entry
                for (m = m_start; m < m_end; ++m) {
                    // (0 - 32), need 32 TLB entry
                    for (k = k_start; k < k_end; ++k) {
#if MATRIX_AB_NEED_TRANSPOSE
                        A_m_k = A[k * M + m];
#else
                        A_m_k = A[m * K + k];
#endif  /* MATRIX_AB_NEED_TRANSPOSE */
                        if (fabs(A_m_k) > (float_t)FLOAT_T_EPSILON) {
                            C_ = &C[m * N];
                            B_ = &B[k * N];
                            // (0 - 512), 512 x 8 = 4096 bytes(1 page size), only need 1 TLB entry,
                            // but in C[m, n], need n_step TLB entry too.
                            for (n = n_start; n < n_end; ++n) {
                                // C[m, n] += A[m, k] * B[k, n];
                                //C[m * N + n] += A_m_k * B[k * N + n];
                                (*C_++) += A_m_k * (*B_++);
                            }
                        }
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
#else
    // 分块循环顺序M, K, N
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
            n_start = 0;
            n_end   = n_start + n_step;
            if (n_end > N)
                n_end = N;
            while (n_start < N) {
                // (0 - 112), need n_step * 2(112 x 2) TLB entry
                for (m = m_start; m < m_end; ++m) {
                    // (0 - 32), need 32 TLB entry
                    for (k = k_start; k < k_end; ++k) {
#if MATRIX_AB_NEED_TRANSPOSE
                        A_m_k = A[k * M + m];
#else
                        A_m_k = A[m * K + k];
#endif  /* MATRIX_AB_NEED_TRANSPOSE */
                        if (fabs(A_m_k) > (float_t)FLOAT_T_EPSILON) {
                            // (0 - 512), 512 x 8 = 4096 bytes(1 page size), only need 1 TLB entry,
                            // but in C[m, n], need n_step TLB entry too.
                            for (n = n_start; n < n_end; ++n) {
                                // C[m, n] += A[m, k] * B[k, n];
                                C[m * N + n] += A_m_k * B[k * N + n];
                            }
                        }
                    }
                }
                n_start += n_step;
                n_end   += n_step;
                if (n_end > N)
                    n_end = N;
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
#endif

#else  /* !CBLAS_ROW_MAJOR */

    float_t B_k_n;
    // fastest (1.27-1.29 seconds)
    m_step = 1024;
    k_step = 64;
    n_step = 64;

#if 0
    m_step = 1024;
    k_step = 1024;
    n_step = 1024;
#endif

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
                // (0 - 32), need 32 TLB entry
                for (k = k_start; k < k_end; ++k) {
                    // (0 - 112), need n_step * 2(112 x 2) TLB entry
                    for (n = n_start; n < n_end; ++n) {
#if MATRIX_AB_NEED_TRANSPOSE
                        B_k_n = B[n + k * N];
#else
                        B_k_n = B[k + n * K];
#endif  /* MATRIX_AB_NEED_TRANSPOSE */
                        if (fabs(B_k_n) > FLOAT_T_EPSILON) {
                            // (0 - 512), 512 x 8 = 4096 bytes(1 page size), only need 1 TLB entry,
                            // but in C[m, n], need n_step TLB entry too.
                            for (m = m_start; m < m_end; ++m) {
                                // C[m, n] += A[m, k] * B[k, n];
                                C[m + n * M] += A[m + k * M] * B_k_n;
                            }
                        }
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
#endif  /* CBLAS_ROW_MAJOR */

}

void matrix_fast_matmult_sse2_4x2(unsigned int M, unsigned int K, unsigned int N,
                                  float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int k_start, k_end;
    unsigned int m_step, n_step, k_step;
    float_t *C_, *B_, *A_;
    float_t A_m_k = (float_t)0.0;

    //
    // matrix multiplication: C = A * B
    //
    // [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
    // [D][E][F] * [I][J]  =  [D*G + E*I + F*K][D*H + E*J + F*L]
    //             [K][L]
    //

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

    // 分块循环顺序K, M, N, 内层循环顺序k, m, n
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
                // k_step (0 - 128),
                // B[k, n] need 128 x 1 TLB entries,
                // A[m, k] need m_step TLB entries only if k_step < nPageSize.
                for (k = k_start; k < k_end; ++k) {
                    // m_step (0 - 4),
                    // C[m, n] need m_step * ((n_step * 8) / 4096) = 4 x 1 TLB entries,
                    // A[m, k] only need m_step x |k_step / nPageSize| = 4 x 1 TLB entries if k_step < nPageSize.
                    for (m = m_start; m < m_end; ++m) {
#if 0
                        A_ = &A[(m + 1) * K + k];
                        _mm_prefetch((char *)A_, _MM_HINT_NTA);

                        A_m_k = A[m * K + k];
                        if (fabs(A_m_k) > (float_t)FLOAT_T_EPSILON) {
                            B_ = &B[k * N + n_start];
                            C_ = &C[m * N + n_start];

                            // n_step (0 - 512),
                            // (nPageSize = 4096 bytes)
                            // B[k, n] need k_step * ((n_step * 8) / nPageSize) = k_step x 1 = 128 x 1 TLB entries,
                            // C[m, n] need m_step * ((n_step * 8) / nPageSize) = m_step x 1 = 4 x 1 TLB entries.
                            for (n = n_start; n < n_end; ++n) {
                                // C[m, n] += A[m, k] * B[k, n];
                                //C[m * N + n] += A_m_k * B[k * N + n];
                                (*C_++) += A_m_k * (*B_++);
                            }
                        }
#elif 1
                        __m128d mm0, mm1, mm2, mm3, rr0, rr1, rr2, rr3;
                        __m128d mma, mmb;
                        A_ = &A[m * K + k];
                        //_mm_prefetch((char *)&A_[k], _MM_HINT_NTA);

                        mma = _mm_load_sd(A_);
                        mmb = _mm_load_sd(A_);
                        mma = _mm_unpacklo_pd(mma, mma);
                        mmb = _mm_unpacklo_pd(mmb, mmb);

                        B_ = &B[k * N + 0];
                        C_ = &C[m * N + 0];

                        for (n = n_start; n < n_end; n += 8) {
                            //_mm_prefetch((char *)&B_[n + 16], _MM_HINT_NTA);
                            //_mm_prefetch((char *)&C_[n + 16], _MM_HINT_NTA);
    #if 1
                            mm0 = _mm_load_pd(&B_[n + 0]);
                            mm1 = _mm_load_pd(&B_[n + 2]);
                            mm2 = _mm_load_pd(&B_[n + 4]);
                            mm3 = _mm_load_pd(&B_[n + 6]);

                            mm0 = _mm_mul_pd(mm0, mma);
                            mm1 = _mm_mul_pd(mm1, mmb);

                            mm2 = _mm_mul_pd(mm2, mma);
                            mm3 = _mm_mul_pd(mm3, mmb);

                            rr0 = _mm_load_pd(&C_[n + 0]);
                            rr1 = _mm_load_pd(&C_[n + 2]);
                            rr2 = _mm_load_pd(&C_[n + 4]);
                            rr3 = _mm_load_pd(&C_[n + 6]);

                            mm0 = _mm_add_pd(mm0, rr0);
                            mm1 = _mm_add_pd(mm1, rr1);
                            mm2 = _mm_add_pd(mm2, rr2);
                            mm3 = _mm_add_pd(mm3, rr3);

                            _mm_store_pd(&C_[n + 0], mm0);
                            _mm_store_pd(&C_[n + 2], mm1);
                            _mm_store_pd(&C_[n + 4], mm2);
                            _mm_store_pd(&C_[n + 6], mm3);
    #else
                            mm0 = _mm_load_pd(&B_[n + 0]);
                            mm1 = _mm_load_pd(&B_[n + 2]);
                            mm2 = _mm_load_pd(&B_[n + 4]);
                            mm3 = _mm_load_pd(&B_[n + 6]);

                            mm0 = _mm_mul_pd(mm0, mma);
                            rr0 = _mm_load_pd(&C_[n + 0]);
                            mm1 = _mm_mul_pd(mm1, mmb);                          
                            rr1 = _mm_load_pd(&C_[n + 2]);

                            mm2 = _mm_mul_pd(mm2, mma);
                            mm0 = _mm_add_pd(mm0, rr0);
                            mm3 = _mm_mul_pd(mm3, mmb);                           
                            mm1 = _mm_add_pd(mm1, rr1);

                            rr2 = _mm_load_pd(&C_[n + 4]);
                            mm2 = _mm_add_pd(mm2, rr2);
                            rr3 = _mm_load_pd(&C_[n + 6]);
                            mm3 = _mm_add_pd(mm3, rr3);

                            _mm_store_pd(&C_[n + 0], mm0);
                            _mm_store_pd(&C_[n + 2], mm1);
                            _mm_store_pd(&C_[n + 4], mm2);
                            _mm_store_pd(&C_[n + 6], mm3);
    #endif
                        }
#elif 0
                        __m128d mm0, mm1, mm2, mm3, rr0, rr1, rr2, rr3;
                        __m128d mma, mmb;
                        A_ = &A[m * K + k];
                        //_mm_prefetch((char *)&A_[k], _MM_HINT_NTA);

                        mma = _mm_load_sd(A_);
                        mmb = _mm_load_sd(A_);
                        mma = _mm_unpacklo_pd(mma, mma);
                        mmb = _mm_unpacklo_pd(mmb, mmb);

                        B_ = &B[k * N + n_start];
                        C_ = &C[m * N + n_start];

                        for (n = (n_end - n_start) / 8; n > 0; B_ += 8, C_ += 8, --n) {
                            //_mm_prefetch((char *)(B_ + 16), _MM_HINT_NTA);
                            //_mm_prefetch((char *)(C_ + 16), _MM_HINT_NTA);

                            mm0 = _mm_load_pd(B_ + 0);
                            mm1 = _mm_load_pd(B_ + 2);
                            mm2 = _mm_load_pd(B_ + 4);
                            mm3 = _mm_load_pd(B_ + 6);

                            mm0 = _mm_mul_pd(mm0, mma);
                            mm1 = _mm_mul_pd(mm1, mmb);

                            mm2 = _mm_mul_pd(mm2, mma);
                            mm3 = _mm_mul_pd(mm3, mmb);

                            rr0 = _mm_load_pd(C_ + 0);
                            rr1 = _mm_load_pd(C_ + 2);
                            rr2 = _mm_load_pd(C_ + 4);
                            rr3 = _mm_load_pd(C_ + 6);

                            mm0 = _mm_add_pd(mm0, rr0);
                            mm1 = _mm_add_pd(mm1, rr1);
                            mm2 = _mm_add_pd(mm2, rr2);
                            mm3 = _mm_add_pd(mm3, rr3);

                            _mm_store_pd(C_ + 0, mm0);
                            _mm_store_pd(C_ + 2, mm1);
                            _mm_store_pd(C_ + 4, mm2);
                            _mm_store_pd(C_ + 6, mm3);
                        }
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
}

#define L1_DCACHE_LINESIZE      64
#define STEP                    (L1_DCACHE_LINESIZE / sizeof(float_t))

#define KERNEL_SEE2_4X2(addr)                                             \
                                                                          \
        __asm { /**-/ prefetchnta byte ptr [eax + (addr) + 24 * 8] /-**/} \
        __asm { movapd      xmm0, xmmword ptr [eax + (addr) + 0 * 8]    } \
        __asm { /**-/ prefetcht0 byte ptr [ecx + (addr) + 24 * 8] /-**/ } \
        __asm { movapd      xmm1, xmmword ptr [eax + (addr) + 2 * 8]    } \
                                                                          \
        __asm { mulpd       xmm0, xmm4                                  } \
        __asm { mulpd       xmm1, xmm5                                  } \
                                                                          \
        __asm { movapd      xmm2, xmmword ptr [eax + (addr) + 4 * 8]    } \
        __asm { movapd      xmm3, xmmword ptr [eax + (addr) + 6 * 8]    } \
                                                                          \
        __asm { addpd       xmm0, xmmword ptr [ecx + (addr) + 0 * 8]    } \
        __asm { addpd       xmm1, xmmword ptr [ecx + (addr) + 2 * 8]    } \
                                                                          \
        __asm { mulpd       xmm2, xmm4                                  } \
        __asm { mulpd       xmm3, xmm5                                  } \
                                                                          \
        __asm { movapd      xmmword ptr [ecx + (addr) + 0 * 8], xmm0    } \
        __asm { movapd      xmmword ptr [ecx + (addr) + 2 * 8], xmm1    } \
                                                                          \
        __asm { addpd       xmm2, xmmword ptr [ecx + (addr) + 4 * 8]    } \
        __asm { addpd       xmm3, xmmword ptr [ecx + (addr) + 6 * 8]    } \
                                                                          \
        __asm { movapd      xmmword ptr [ecx + (addr) + 4 * 8], xmm2    } \
        __asm { movapd      xmmword ptr [ecx + (addr) + 6 * 8], xmm3    }
                                                            
void matmult_s_tiling_sse2_4x2(unsigned int M, unsigned int K, unsigned int N,
                               float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m2, n2, k2;
    /*
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int k_start, k_end;
    //*/
    unsigned int m_step, n_step, k_step;
    float_t *mul1, *mul2, *res;
    float_t *RESTRICT rres;
    float_t *RESTRICT rmul1;
    float_t *RESTRICT rmul2;

    mul1 = A;
    mul2 = B;
    res  = C;

    m_step = 4;
    k_step = 128;
    n_step = 512;

    m_step = 4;
    k_step = 128;
    n_step = 512;

#if 0
    m_step = 4;
    k_step = 256;
    n_step = 512;
#endif

#if 0
    m_step = 2;
    k_step = 128;
    n_step = 1024;
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

    //__asm { emms }

    for (k = 0; k < K; k += k_step) {
        for (m = 0; m < M; m += m_step) {
          for (n = 0; n < N; n += n_step) {
            //for (m = 0; m < M; m += m_step) {
                for (k2 = 0, rmul2 = &mul2[k * N + n]; k2 < k_step; ++k2, rmul2 += N) {
                    //_mm_prefetch((char *)&rmul1[8], _MM_HINT_NTA);
                    for (m2 = 0, rres = &res[m * N + n], rmul1 = &mul1[m * K + k]; m2 < m_step; ++m2, rres += N, rmul1 += K) {
#if 0
                        __m128d mm0, mm1, rr0, rr1;
                        __m128d mma, mmb;
                        //_mm_prefetch((char *)&rmul1[K], _MM_HINT_NTA);
                        _mm_prefetch((char *)&rmul1[K], _MM_HINT_T0);
                        mma = _mm_load_sd(&rmul1[k2]);
                        mmb = _mm_load_sd(&rmul1[k2]);
                        mma = _mm_unpacklo_pd(mma, mma);
                        mmb = _mm_unpacklo_pd(mmb, mmb);
                        for (n2 = 0; n2 < n_step; n2 += 8) {
                            //_mm_prefetch((char *)&rmul2[n2 + 16], _MM_HINT_NTA);
                            //_mm_prefetch((char *)&rres[n2 + 16],  _MM_HINT_NTA);

                            mm0 = _mm_load_pd(&rmul2[n2 + 0]);
                            mm1 = _mm_load_pd(&rmul2[n2 + 2]);

                            mm0 = _mm_mul_pd(mm0, mma);
                            mm1 = _mm_mul_pd(mm1, mmb);

                            rr0 = _mm_load_pd(&rres[n2 + 0]);
                            rr1 = _mm_load_pd(&rres[n2 + 2]);

                            rr0 = _mm_add_pd(rr0, mm0);
                            rr1 = _mm_add_pd(rr1, mm1);

                            _mm_store_pd(&rres[n2 + 0], rr0);
                            _mm_store_pd(&rres[n2 + 2], rr1);

                            mm0 = _mm_load_pd(&rmul2[n2 + 4]);
                            mm1 = _mm_load_pd(&rmul2[n2 + 6]);

                            mm0 = _mm_mul_pd(mm0, mma);
                            mm1 = _mm_mul_pd(mm1, mmb);

                            rr0 = _mm_load_pd(&rres[n2 + 4]);
                            rr1 = _mm_load_pd(&rres[n2 + 6]);

                            rr0 = _mm_add_pd(rr0, mm0);
                            rr1 = _mm_add_pd(rr1, mm1);

                            _mm_store_pd(&rres[n2 + 4], rr0);
                            _mm_store_pd(&rres[n2 + 6], rr1);
                        }
#elif 0
                        __m128d mm0, mm1, mm2, mm3, rr0, rr1, rr2, rr3;
                        __m128d mma, mmb;
                        _mm_prefetch((char *)&rmul1[K], _MM_HINT_NTA);
                        mma = _mm_load_sd(&rmul1[k2]);
                        mmb = _mm_load_sd(&rmul1[k2]);
                        mma = _mm_unpacklo_pd(mma, mma);
                        mmb = _mm_unpacklo_pd(mmb, mmb);
                        for (n2 = 0; n2 < n_step; n2 += 8) {
                            //_mm_prefetch((char *)&rmul2[n2 + 16], _MM_HINT_NTA);
                            //_mm_prefetch((char *)&rres[n2 + 16],  _MM_HINT_NTA);

                            mm0 = _mm_load_pd(&rmul2[n2 + 0]);
                            mm1 = _mm_load_pd(&rmul2[n2 + 2]);
                            mm2 = _mm_load_pd(&rmul2[n2 + 4]);
                            mm3 = _mm_load_pd(&rmul2[n2 + 6]);

                            mm0 = _mm_mul_pd(mm0, mma);
                            mm1 = _mm_mul_pd(mm1, mmb);

                            mm2 = _mm_mul_pd(mm2, mma);
                            mm3 = _mm_mul_pd(mm3, mmb);

                            rr0 = _mm_load_pd(&rres[n2 + 0]);
                            rr1 = _mm_load_pd(&rres[n2 + 2]);
                            rr2 = _mm_load_pd(&rres[n2 + 4]);
                            rr3 = _mm_load_pd(&rres[n2 + 6]);

                            mm0 = _mm_add_pd(mm0, rr0);
                            mm1 = _mm_add_pd(mm1, rr1);
                            mm2 = _mm_add_pd(mm2, rr2);
                            mm3 = _mm_add_pd(mm3, rr3);

                            _mm_store_pd(&rres[n2 + 0], mm0);
                            _mm_store_pd(&rres[n2 + 2], mm1);
                            _mm_store_pd(&rres[n2 + 4], mm2);
                            _mm_store_pd(&rres[n2 + 6], mm3);
                        }
#else
                        float_t *pmul1, *pmul2, *pres;
                        pmul1 = &rmul1[k2];
                        pmul2 = &rmul2[0];
                        pres  = &rres[0];
                        n2 = 0;
                        __asm {
                            push        edi;
                            push        eax;
                            push        ecx;
                            push        ebx;

                            mov         edi, pmul1;
                            mov         eax, pmul2;
                            mov         ecx, pres;
                            //mov         ebx, K
                            movsd       xmm4, qword ptr [edi];
                            movapd      xmm5, xmm4;
                            //prefetcht0  byte ptr [edi + ebx];
                            mov         ebx, n_step;
                            unpcklpd    xmm4, xmm4;
                            unpcklpd    xmm5, xmm5;

                            ALIGN_16
L01:
                            // for (n2 = 0; n2 < n_step; n2 += 8) {

                            #define SSE2_PREFETCH_SIZE  ((15 * 8 + 4) * FLOAT_SIZE)

    #if 1
                            ///////////////////////////////////////////////

                            //prefetcht0  byte ptr [eax + SSE2_PREFETCH_SIZE + 0 * 8];
                            //prefetcht0  byte ptr [ecx + SSE2_PREFETCH_SIZE + 0 * 8];

                            //movq        mm2, qword ptr [eax + 32 * 8];
                            movapd      xmm0, xmmword ptr [eax + 0 * 8];
                            movapd      xmm1, xmmword ptr [eax + 2 * 8];

                            mulpd       xmm0, xmm4;
                            movapd      xmm2, xmmword ptr [eax + 4 * 8];
                            mulpd       xmm1, xmm5;
                            movapd      xmm3, xmmword ptr [eax + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + 0 * 8];
                            mulpd       xmm2, xmm4;
                            //prefetchnta byte ptr [eax + SSE2_PREFETCH_SIZE + 8 * 8];
                            //prefetcht0  byte ptr [eax + SSE2_PREFETCH_SIZE + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + 2 * 8];
                            mulpd       xmm3, xmm5;
                            //prefetcht0  byte ptr [ecx + SSE2_PREFETCH_SIZE + 0 * 8];

                            addpd       xmm2, xmmword ptr [ecx + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + 6 * 8];

                            ///*
                            movapd      xmmword ptr [ecx + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + 2 * 8], xmm1;
                            movapd      xmmword ptr [ecx + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + 6 * 8], xmm3;
                            //*/

                            /*
                            movdqa      xmmword ptr [ecx + 0 * 8], xmm0;
                            movdqa      xmmword ptr [ecx + 2 * 8], xmm1;
                            movdqa      xmmword ptr [ecx + 4 * 8], xmm2;
                            movdqa      xmmword ptr [ecx + 6 * 8], xmm3;
                            //*/

                            ///////////////////////////////////////////////

                            //prefetcht0  byte ptr [eax + SSE2_PREFETCH_SIZE +  8 * 8];
                            //prefetcht0  byte ptr [ecx + SSE2_PREFETCH_SIZE +  8 * 8];

                            //movq        mm3, qword ptr [eax + 40 * 8];
                            movapd      xmm0, xmmword ptr [eax +  8 * 8];
                            movapd      xmm1, xmmword ptr [eax + 10 * 8];

                            mulpd       xmm0, xmm4;
                            movapd      xmm2, xmmword ptr [eax + 12 * 8];
                            mulpd       xmm1, xmm5;
                            movapd      xmm3, xmmword ptr [eax + 14 * 8];

                            addpd       xmm0, xmmword ptr [ecx +  8 * 8];
                            mulpd       xmm2, xmm4;
                            //prefetchnta byte ptr [eax + SSE2_PREFETCH_SIZE + 8 * 8];
                            //prefetcht0  byte ptr [eax + SSE2_PREFETCH_SIZE + 8 * 8];
                            addpd       xmm1, xmmword ptr [ecx + 10 * 8];
                            mulpd       xmm3, xmm5;
                            //prefetcht0  byte ptr [ecx + SSE2_PREFETCH_SIZE + 8 * 8];

                            addpd       xmm2, xmmword ptr [ecx + 12 * 8];
                            addpd       xmm3, xmmword ptr [ecx + 14 * 8];

                            ///*
                            movapd      xmmword ptr [ecx +  8 * 8], xmm0;
                            movapd      xmmword ptr [ecx + 10 * 8], xmm1;
                            movapd      xmmword ptr [ecx + 12 * 8], xmm2;
                            movapd      xmmword ptr [ecx + 14 * 8], xmm3;
                            //*/

                            /*
                            movdqa      xmmword ptr [ecx +  8 * 8], xmm0;
                            movdqa      xmmword ptr [ecx + 10 * 8], xmm1;
                            movdqa      xmmword ptr [ecx + 12 * 8], xmm2;
                            movdqa      xmmword ptr [ecx + 14 * 8], xmm3;
                            //*/

                            //prefetcht0  byte ptr [ecx + SSE2_PREFETCH_SIZE + 8 * 8];

                            ///////////////////////////////////////////////
    #else
                            ///////////////////////////////////////////////

        #if 1
                            ///*
                            KERNEL_SEE2_4X2( 0 * 8);
                            KERNEL_SEE2_4X2( 8 * 8);
                            KERNEL_SEE2_4X2(16 * 8);
                            KERNEL_SEE2_4X2(24 * 8);
                            KERNEL_SEE2_4X2(32 * 8);
                            KERNEL_SEE2_4X2(40 * 8);
                            KERNEL_SEE2_4X2(48 * 8);
                            KERNEL_SEE2_4X2(56 * 8);
                            //*/
        #else
                            #undef  BASE_ADDR
                            #define BASE_ADDR (0 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

                            #undef  BASE_ADDR
                            #define BASE_ADDR (8 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

                            #undef  BASE_ADDR
                            #define BASE_ADDR (16 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

                            #undef  BASE_ADDR
                            #define BASE_ADDR (24 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

            #if 0
                            #undef  BASE_ADDR
                            #define BASE_ADDR (32 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

                            #undef  BASE_ADDR
                            #define BASE_ADDR (40 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

                            #undef  BASE_ADDR
                            #define BASE_ADDR (48 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;

                            #undef  BASE_ADDR
                            #define BASE_ADDR (56 * 8)

                            //prefetchnta byte ptr [eax + (BASE_ADDR) + 16 * 8];
                            movapd      xmm0, xmmword ptr [eax + (BASE_ADDR) + 0 * 8];
                            //prefetchnta byte ptr [ecx + (BASE_ADDR) + 16 * 8];
                            movapd      xmm1, xmmword ptr [eax + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm0, xmm4;
                            mulpd       xmm1, xmm5;

                            movapd      xmm2, xmmword ptr [eax + (BASE_ADDR) + 4 * 8];
                            movapd      xmm3, xmmword ptr [eax + (BASE_ADDR) + 6 * 8];

                            addpd       xmm0, xmmword ptr [ecx + (BASE_ADDR) + 0 * 8];
                            addpd       xmm1, xmmword ptr [ecx + (BASE_ADDR) + 2 * 8];

                            mulpd       xmm2, xmm4;
                            mulpd       xmm3, xmm5;

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 0 * 8], xmm0;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 2 * 8], xmm1;

                            addpd       xmm2, xmmword ptr [ecx + (BASE_ADDR) + 4 * 8];
                            addpd       xmm3, xmmword ptr [ecx + (BASE_ADDR) + 6 * 8];

                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 4 * 8], xmm2;
                            movapd      xmmword ptr [ecx + (BASE_ADDR) + 6 * 8], xmm3;
            #endif
                            ///////////////////////////////////////////////
        #endif

    #endif

    #if 0
                            add         eax, 64;
                            add         ecx, 64;
                            sub         ebx, 8;
    #elif 1
                            add         eax, 64 * 2;
                            add         ecx, 64 * 2;
                            sub         ebx,  8 * 2;
    #elif 0
                            add         eax, 64 * 4;
                            add         ecx, 64 * 4;
                            sub         ebx,  8 * 4;
    #else
                            add         eax, 64 * 8;
                            add         ecx, 64 * 8;
                            sub         ebx,  8 * 8;
    #endif
                            jnz         L01
                            // }
                            pop         ebx;
                            pop         ecx;
                            pop         eax;
                            pop         edi;
                        }
#endif
                    }
                }
#if 0
                __asm { sfence }
#endif
            }
        }
    }

    //__asm { sfence }
    //__asm { emms }
}

void matmult_s_tiling_sse2(unsigned int M, unsigned int K, unsigned int N,
                           float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m2, n2, k2;
    /*
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int k_start, k_end;
    //*/
    unsigned int m_step, n_step, k_step;
    float_t *mul1, *mul2, *res;
    float_t *RESTRICT rres;
    float_t *RESTRICT rmul1;
    float_t *RESTRICT rmul2;

    mul1 = A;
    mul2 = B;
    res  = C;

    m_step = 4;
    k_step = 128;
    n_step = 512;

    if (m_step > M)
        m_step = M;
    if (k_step > K)
        k_step = K;
    if (n_step > N)
        n_step = N;

    for (k = 0; k < K; k += k_step) {
        for (m = 0; m < M; m += m_step) {
            for (n = 0; n < N; n += n_step) {
                for (k2 = 0, rmul2 = &mul2[k * N + n]; k2 < k_step; ++k2, rmul2 += N) {
                    //_mm_prefetch((char *)&rmul1[8], _MM_HINT_NTA);
                    for (m2 = 0, rres = &res[m * N + n], rmul1 = &mul1[m * K + k]; m2 < m_step; ++m2, rres += N, rmul1 += K) {
                        __m128d m2, r2;
                        __m128d m1d;
                        _mm_prefetch((char *)&rmul1[K], _MM_HINT_NTA);
                        m1d = _mm_load_sd(&rmul1[k2]);
                        m1d = _mm_unpacklo_pd(m1d, m1d);
                        for (n2 = 0; n2 < n_step; n2 += 8) {
                            //_mm_prefetch((char *)&rmul2[n2 + 16], _MM_HINT_NTA);
                            //_mm_prefetch((char *)&rres[n2 + 16],  _MM_HINT_NTA);

                            m2 = _mm_load_pd(&rmul2[n2]);
                            r2 = _mm_load_pd(&rres[n2]);
                            _mm_store_pd(&rres[n2], _mm_add_pd(_mm_mul_pd(m2, m1d), r2));

                            m2 = _mm_load_pd(&rmul2[n2 + 2]);
                            r2 = _mm_load_pd(&rres[n2 + 2]);
                            _mm_store_pd(&rres[n2 + 2], _mm_add_pd(_mm_mul_pd(m2, m1d), r2));

                            m2 = _mm_load_pd(&rmul2[n2 + 4]);
                            r2 = _mm_load_pd(&rres[n2 + 4]);
                            _mm_store_pd(&rres[n2 + 4], _mm_add_pd(_mm_mul_pd(m2, m1d), r2));

                            m2 = _mm_load_pd(&rmul2[n2 + 6]);
                            r2 = _mm_load_pd(&rres[n2 + 6]);
                            _mm_store_pd(&rres[n2 + 6], _mm_add_pd(_mm_mul_pd(m2, m1d), r2));
                        }
                    }
                }
            }
        }
    }
}

void matmult_s_tiling_sse2_1(unsigned int M, unsigned int K, unsigned int N,
                             float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m2, n2, k2;
    unsigned int m_step, n_step, k_step;
    float_t *mul1, *mul2, *res;
    float_t *RESTRICT rres;
    float_t *RESTRICT rmul1;
    float_t *RESTRICT rmul2;

    mul1 = A;
    mul2 = B;
    res  = C;

    m_step = 8;
    k_step = 64;
    n_step = 512;

    if (m_step > M)
        m_step = M;
    if (k_step > K)
        k_step = K;
    if (n_step > N)
        n_step = N;

    for (m = 0; m < M; m += m_step) {
        for (n = 0; n < N; n += n_step) {
            for (k = 0; k < K; k += k_step) {
                for (m2 = 0, rres = &res[m * N + n], rmul1 = &mul1[m * K + k]; m2 < m_step; ++m2, rres += N, rmul1 += K) {
                    _mm_prefetch((char *)&rmul1[8], _MM_HINT_NTA);
                    for (k2 = 0, rmul2 = &mul2[k * N + n]; k2 < k_step; ++k2, rmul2 += N) {
                        __m128d m1d = _mm_load_sd(&rmul1[k2]);
                        m1d = _mm_unpacklo_pd(m1d, m1d);
                        for (n2 = 0; n2 < n_step; n2 += 2) {
                            if ((n2 & 0x00000007UL) == 0) {
                                _mm_prefetch((char *)&rmul2[n2 + 8], _MM_HINT_NTA);
                                _mm_prefetch((char *)&rres[n2 + 8],  _MM_HINT_NTA);
                            }
                            __m128d m2 = _mm_load_pd(&rmul2[n2]);
                            __m128d r2 = _mm_load_pd(&rres[n2]);
                            _mm_store_pd(&rres[n2], _mm_add_pd(_mm_mul_pd(m2, m1d), r2));
                        }
                    }
                }
            }
        }
    }
}

void matmult_s_tiling_sse2_0(unsigned int M, unsigned int K, unsigned int N,
                             float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m2, n2, k2;
    //unsigned int m_step, n_step, k_step;
    float_t *mul1, *mul2, *res;
    float_t *RESTRICT rres;
    float_t *RESTRICT rmul1;
    float_t *RESTRICT rmul2;

    mul1 = A;
    mul2 = B;
    res  = C;

    for (m = 0; m < M; m += STEP) {
        for (n = 0; n < N; n += STEP) {
            for (k = 0; k < K; k += STEP) {
                for (m2 = 0, rres = &res[m * N + n], rmul1 = &mul1[m * K + k]; m2 < STEP; ++m2, rres += N, rmul1 += K) {
                    _mm_prefetch((char *)&rmul1[8], _MM_HINT_NTA);
                    for (k2 = 0, rmul2 = &mul2[k * N + n]; k2 < STEP; ++k2, rmul2 += N) {
                        __m128d m1d = _mm_load_sd(&rmul1[k2]);
                        m1d = _mm_unpacklo_pd(m1d, m1d);
                        for (n2 = 0; n2 < STEP; n2 += 2) {
                            __m128d m2 = _mm_load_pd(&rmul2[n2]);
                            __m128d r2 = _mm_load_pd(&rres[n2]);
                            _mm_store_pd(&rres[n2], _mm_add_pd(_mm_mul_pd(m2, m1d), r2));
                        }
                    }
                }
            }
        }
    }
}

void serial_matmult_sse(unsigned int M, unsigned int K, unsigned int N,
                        float_t *A, float_t *B, float_t *C)
{
#if 1
    float_t temp;
    unsigned int m, k, n;
    register float_t *A_ptr;
    register float_t *C_ptr, *B_ptr;
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            // C[m][n] = 0.0;
            C[m * N + n] = (float_t)0.0;
        }
        A_ptr = (float_t *)&(A[m * K]);
        for (k = 0; k < K; ++k) {
            temp = *(A_ptr + k);
            //if (temp > (float_t)FLOAT_T_EPSILON || temp < (float_t)-FLOAT_T_EPSILON) {
            if (fabs(temp) > (float_t)FLOAT_T_EPSILON) {
                C_ptr = (float_t *)&(C[m * N]);
                B_ptr = (float_t *)&(B[k * N]);

#if (_MULT_SSE2_MODE == 1)
                for (n = 0; n < N; ++n) {
                    // C[m * N + n] += temp * B[k * N + n];
                    (*C_ptr++) += temp * (*B_ptr++);
                }

#elif (_MULT_SSE2_MODE == 2)
                __m128d tmp;
                __m128d mm1, mm2, mm3, mm4;
#if USE_SEE2_WRITE_PREFRTCH
                __m128d mma, mmb;
#endif
#if (DOUBLE_NUMS_PER_LOOP == 16)
                __m128d mm5, mm6, mm7, mm8;
#endif
                tmp = _mm_setzero_pd();
                tmp = _mm_loadl_pd(tmp, &temp);
                tmp = _mm_shuffle_pd(tmp, tmp, 0);

                for (n = N / DOUBLE_NUMS_PER_LOOP;
                     n > 0;
                     C_ptr += DOUBLE_NUMS_PER_LOOP, B_ptr += DOUBLE_NUMS_PER_LOOP, --n) {
#if USE_SEE2_WRITE_PREFRTCH
                    // prefetchnta [addr + offset]
                    _mm_prefetch(((char *)B_ptr + MM_PREFETCH_OFFSET), _MM_HINT_NTA);
                    _mm_prefetch(((char *)C_ptr + MM_PREFETCH_OFFSET), _MM_HINT_NTA);

                    mm1 = _mm_load_pd(B_ptr   );
                    mm2 = _mm_load_pd(B_ptr + 2);

                    mma = _mm_load_pd(C_ptr);
                    mmb = _mm_load_pd(C_ptr + 2);

                    mm1 = _mm_mul_pd(mm1, tmp);
                    mm2 = _mm_mul_pd(mm2, tmp);

                    mm3 = _mm_load_pd(B_ptr + 4);
                    mm4 = _mm_load_pd(B_ptr + 6);

                    mm3 = _mm_mul_pd(mm3, tmp);
                    mm4 = _mm_mul_pd(mm4, tmp);

                    mm1 = _mm_add_pd(mm1, mma);
                    mm2 = _mm_add_pd(mm2, mmb);

                    mm3 = _mm_add_pd(mm3, *((__m128d *)(C_ptr + 4)));
                    mm4 = _mm_add_pd(mm4, *((__m128d *)(C_ptr + 6)));

                    _mm_store_pd(C_ptr,     mm1);
                    _mm_store_pd(C_ptr + 2, mm2);

#if (DOUBLE_NUMS_PER_LOOP == 16)
                    _mm_prefetch(((char *)B_ptr + MM_PREFETCH_OFFSET + 64), _MM_HINT_NTA);
                    _mm_prefetch(((char *)C_ptr + MM_PREFETCH_OFFSET + 64), _MM_HINT_NTA);

                    mm5 = _mm_load_pd(B_ptr + 8);
                    mm6 = _mm_load_pd(B_ptr + 10);
#endif
                    _mm_store_pd(C_ptr + 4, mm3);
                    _mm_store_pd(C_ptr + 6, mm4);

#if (DOUBLE_NUMS_PER_LOOP == 16)
                    mma = _mm_load_pd(C_ptr + 8);
                    mmb = _mm_load_pd(C_ptr + 10);

                    mm5 = _mm_mul_pd(mm5, tmp);
                    mm6 = _mm_mul_pd(mm6, tmp);

                    mm7 = _mm_load_pd(B_ptr + 12);
                    mm8 = _mm_load_pd(B_ptr + 14);

                    mm5 = _mm_add_pd(mm5, mma);
                    mm6 = _mm_add_pd(mm6, mmb);

                    mm7 = _mm_mul_pd(mm7, tmp);
                    mm8 = _mm_mul_pd(mm8, tmp);

                    _mm_store_pd(C_ptr + 8,  mm5);
                    _mm_store_pd(C_ptr + 10, mm6);

                    mm7 = _mm_add_pd(mm7, *((__m128d *)(C_ptr + 12)));
                    mm8 = _mm_add_pd(mm8, *((__m128d *)(C_ptr + 14)));

                    _mm_store_pd(C_ptr + 12, mm7);
                    _mm_store_pd(C_ptr + 14, mm8);
#endif
#else /* !USE_SEE2_WRITE_PREFRTCH */
                    // prefetchnta [addr + offset]
                    //_mm_prefetch(((char *)B_ptr + MM_PREFETCH_OFFSET), _MM_HINT_NTA);
                    //_mm_prefetch(((char *)C_ptr + MM_PREFETCH_OFFSET), _MM_HINT_NTA);

                    mm1 = _mm_load_pd(B_ptr   );
                    mm2 = _mm_load_pd(B_ptr + 2);

                    mm3 = _mm_load_pd(B_ptr + 4);
                    mm4 = _mm_load_pd(B_ptr + 6);

                    mm1 = _mm_mul_pd(mm1, tmp);
                    mm2 = _mm_mul_pd(mm2, tmp);

                    mm3 = _mm_mul_pd(mm3, tmp);
                    mm4 = _mm_mul_pd(mm4, tmp);

                    mm1 = _mm_add_pd(mm1, *((__m128d *)(C_ptr))   );
                    mm2 = _mm_add_pd(mm2, *((__m128d *)(C_ptr + 2)));

                    mm3 = _mm_add_pd(mm3, *((__m128d *)(C_ptr + 4)));
                    mm4 = _mm_add_pd(mm4, *((__m128d *)(C_ptr + 6)));

                    _mm_store_pd(C_ptr,     mm1);
                    _mm_store_pd(C_ptr + 2, mm2);

#if (DOUBLE_NUMS_PER_LOOP == 16)
                    //_mm_prefetch(((char *)B_ptr + MM_PREFETCH_OFFSET + 64), _MM_HINT_NTA);
                    //_mm_prefetch(((char *)C_ptr + MM_PREFETCH_OFFSET + 64), _MM_HINT_NTA);

                    mm5 = _mm_load_pd(B_ptr + 8);
                    mm6 = _mm_load_pd(B_ptr + 10);
#endif
                    _mm_store_pd(C_ptr + 4, mm3);
                    _mm_store_pd(C_ptr + 6, mm4);

#if (DOUBLE_NUMS_PER_LOOP == 16)
                    mm5 = _mm_mul_pd(mm5, tmp);
                    mm6 = _mm_mul_pd(mm6, tmp);

                    mm7 = _mm_load_pd(B_ptr + 12);
                    mm8 = _mm_load_pd(B_ptr + 14);

                    mm5 = _mm_add_pd(mm5, *((__m128d *)(C_ptr +  8)));
                    mm6 = _mm_add_pd(mm6, *((__m128d *)(C_ptr + 10)));

                    mm7 = _mm_mul_pd(mm7, tmp);
                    mm8 = _mm_mul_pd(mm8, tmp);

                    _mm_store_pd(C_ptr + 8,  mm5);
                    _mm_store_pd(C_ptr + 10, mm6);

                    mm7 = _mm_add_pd(mm7, *((__m128d *)(C_ptr + 12)));
                    mm8 = _mm_add_pd(mm8, *((__m128d *)(C_ptr + 14)));

                    _mm_store_pd(C_ptr + 12, mm7);
                    _mm_store_pd(C_ptr + 14, mm8);
#endif
        #endif  /* !USE_SEE2_WRITE_PREFRTCH */
                }
    #else
                for (n = 0; n < N; ++n) {
                    (*C_ptr++) += temp * (*B_ptr++);
                }
    #endif  /* _MULT_SSE2_MODE other value */
            }
        }
    }
#else
    double temp;
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n)
            C[m * N + n] = 0.0;
        for (k = 0; k < K; ++k) {
            temp = A[m * K + k];
            //if (temp > (float_t)DOUBLE_EPSINON || temp < (float_t)-DOUBLE_EPSINON) {
            if (fabs(temp) > (float_t)DOUBLE_EPSINON) {
                for (n = 0; n < N; ++n)
                    C[m * N + n] += temp * B[k * N + n];
            }
        }
    }
#endif
}

void serial_matmult(unsigned int M, unsigned int K, unsigned int N,
                    float_t *A, float_t *B, float_t *C)
{
    float_t tempA;
    unsigned int m, k, n;
    register float_t *A_ptr;
    register float_t *C_ptr, *B_ptr;
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            // C[m][n] = 0.0;
            C[m * N + n] = (float_t)0.0;
        }
        // A_ptr = &A[m][0];
        A_ptr = (float_t *)&(A[m * K + 0]);
        for (k = 0; k < K; ++k) {
            // tempA = A[m][k];
            tempA = *(A_ptr + k);
            //tempA = A[m * K + k];
            //if (tempA > (float_t)FLOAT_T_EPSILON || tempA < (float_t)-FLOAT_T_EPSILON) {
            if (fabs(tempA) > (float_t)FLOAT_T_EPSILON) {
                C_ptr = (float_t *)&(C[m * N + 0]);
                B_ptr = (float_t *)&(B[k * N + 0]);

                for (n = 0; n < N; ++n) {
                    // C[m * N + n] += tempA * B[k * N + n];
                    (*C_ptr++) += tempA * (*B_ptr++);
                    //C[m * N + n] += tempA * B[k * N + n];
                }
            }
        }
    }
}

void matmult_s_row_tiling_MxN_K_transB(unsigned int M, unsigned int K, unsigned int N,
                                       float_t *A, float_t *B, float_t *C)
{
    unsigned int m, n, k;
    unsigned int m_start, m_end;
    unsigned int n_start, n_end;
    unsigned int k_start, k_end;
    unsigned int m_step, n_step, k_step;
    float_t *C_ = NULL, *B_, *A_;
    float_t C_m_n;

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
                        if (k_step >= K)
                            C_m_n = (float_t)0.0;
                        else
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

void matrix_matmult_test(unsigned int M, unsigned int K, unsigned int N)
{
    const unsigned int kAlignment = 128;
    stop_watch stopWatch;
    double elapsedTime1, elapsedTime2, elapsedTime3, elapsedTime4;
    int err_nums = 0, diff_nums = 0;
    bool verify_ok = false;

    float_t *A, *B, *C1, *C2, *C3, *C4;
    float_t *tmp1, *tmp2;;
    float_t alpha = 1.0, beta = 0.0;

    printf("matrix_matmult_test().\n\n");
    printf("M = %d, K = %d, N = %d\n\n", M, K, N);

#if USE_LARGE_PAGES
    const int addr_offset = 128;
    const unsigned int total_matrix_cnt = 6 - 3;
    int total_offset_cnt;
    if (total_matrix_cnt > 1)
        total_offset_cnt = total_matrix_cnt * (total_matrix_cnt - 1) / 2 + 2;
    else
        total_offset_cnt = 2;
    size_t alloc_size = (M * K + K * N + M * N * 4) * sizeof(float_t) + addr_offset * total_offset_cnt;
    void *map_address = NULL;
    int matrix_idx = total_matrix_cnt - 1;
    char *tmp_address;
    tmp1 = tmp2 = NULL;
    map_address = huge_tlb_malloc(alloc_size, &alloc_size);
    if (map_address) {
        printf("map_address = 0x%08X, alloc_size = 0x%08X (%d) byte(s).\n\n", map_address, alloc_size, alloc_size);
        tmp_address = (char *)map_address;
        A  = (float_t *)tmp_address + addr_offset * matrix_idx;
        matrix_idx--;
        if (matrix_idx < 0) matrix_idx = 0;

        tmp_address += M * K * sizeof(float_t) + addr_offset * matrix_idx;
        B  = (float_t *)tmp_address;
        matrix_idx--;
        if (matrix_idx < 0) matrix_idx = 0;

        tmp_address += K * N * sizeof(float_t) + addr_offset * matrix_idx;
        C1 = (float_t *)tmp_address;
        //matrix_idx--;

        tmp_address += M * N * sizeof(float_t) + addr_offset * matrix_idx;
        C2 = (float_t *)tmp_address;
        //matrix_idx--;

        tmp_address += M * N * sizeof(float_t) + addr_offset * matrix_idx;
        C3 = (float_t *)tmp_address;
        //matrix_idx--;

        tmp_address += M * N * sizeof(float_t) + addr_offset * matrix_idx;
        C4 = (float_t *)tmp_address;
        //matrix_idx--;

        tmp_address += M * N * sizeof(float_t) + addr_offset * matrix_idx;

        ::fflush(stdout);

        //Sleep(1000);
        //goto MATRIX_MULT_EXIT;
    }
    else
        return;
#else
    // matrix A
    A  = matrix_malloc(M, K, kAlignment);
    if (A == NULL) {
        // malloc matrix A failure
        printf("matrix_malloc() matrix A failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    tmp1 = matrix_malloc(16, 8, 16);
    if (tmp1 == NULL) {
        goto MATRIX_MULT_EXIT;
    }

    // matrix B
    B  = matrix_malloc(K, N, kAlignment);
    if (B == NULL) {
        // malloc matrix B failure
        printf("matrix_malloc() matrix B failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    tmp2 = matrix_malloc(16, 8, 16);
    if (tmp2 == NULL) {
        goto MATRIX_MULT_EXIT;
    }

    // matrix C1
    C1 = matrix_malloc(M, N, kAlignment);
    if (C1 == NULL) {
        // malloc matrix C1 failure
        printf("matrix_malloc() matrix C1 failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C2
    C2 = matrix_malloc(M, N, kAlignment);
    if (C2 == NULL) {
        // malloc matrix C2 failure
        printf("matrix_malloc() matrix C2 failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C3
    C3 = matrix_malloc(M, N, kAlignment);
    if (C3 == NULL) {
        // malloc matrix C failure
        printf("matrix_malloc() matrix C3 failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C4
    C4 = matrix_malloc(M, N, kAlignment);
    if (C4 == NULL) {
        // malloc matrix C failure
        printf("matrix_malloc() matrix C4 failure.\n");
        goto MATRIX_MULT_EXIT;
    }
#endif  /* USE_LARGE_PAGES */

    // set fixed random seed
    // (设置固定的随机种子, 是为了保证矩阵的初始化数据是一致的, 这样不同的库比较的时候比较公平)
    ::srand(RAND_FIXED_SEED);

    // init elements for all matrixs
    matrix_init_elements(A,  M, K, MatInitRands);
    matrix_init_elements(B,  K, N, MatInitRands);
    matrix_init_elements(C1, M, N, MatInitZeros);
    matrix_init_elements(C2, M, N, MatInitZeros);
    matrix_init_elements(C3, M, N, MatInitZeros);
    matrix_init_elements(C4, M, N, MatInitZeros);

    elapsedTime1 = 0.0;

    /**********************************
     *     matrix_fast_matmult()      *
     **********************************/

    stopWatch.start();
    matrix_fast_matmult(M, K, N, A, B, C1);
    stopWatch.stop();
    elapsedTime2 = stopWatch.getMillisec();

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matrix_fast_matmult: Step1], elapsed time: %0.2f ms\n\n", elapsedTime1);
    printf("[matrix_fast_matmult: Step2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#else
    printf("[matrix_fast_matmult], elapsed time: %0.2f ms\n\n", elapsedTime2);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if 0
    /********************************************************
     *     matmult_s_row_tiling_N_sse2_4x1()      *
     ********************************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_N_sse2_4x1(M, N, K, &alpha, A, K, B, N, &beta, C3, N);
    stopWatch.stop();
    elapsedTime1 = stopWatch.getMillisec();

    printf("[matmult_s_row_tiling_N_sse2_4x1], elapsed time: %0.2f ms\n\n", elapsedTime1);
#endif

    /*************************************
     *     gemm_kernel_2x4_penryn()      *
     *************************************/
    float_t *C3b = C3;
    matrix_init_elements(C3b, M, N, MatInitZeros);

    // 先转置矩阵A
    //matrix_fast_transpose_NxN((float_t *)A, M, K);
    matrix_fast_transpose_NxN((float_t *)B, K, N);

    stopWatch.start();
    //gemm_kernel_2x4_penryn(M, N, K, alpha, A, K, B, N, beta, C3b, N, 0);
    gemm_kernel_2x4_penryn(M, N, K, alpha, A, K, B, K, beta, C3b, M, 0);
    stopWatch.stop();
    elapsedTime1 = stopWatch.getMillisec();

    // 计算完以后再转置(还原)矩阵A
    //matrix_fast_transpose_NxN((float_t *)A, K, M);
    matrix_fast_transpose_NxN((float_t *)B, N, K);

    matrix_fast_transpose_NxN((float_t *)C3, N, M);

    //C3 -= M * N;
    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

    printf("[gemm_kernel_2x4_penryn], elapsed time: %0.2f ms\n\n", elapsedTime1);

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*
    printf("C1[0] = %0.5f, C3[0] = %0.5f\n", C1[0], C3[0]);
    printf("C1[1] = %0.5f, C3[1] = %0.5f\n", C1[1], C3[1]);
    printf("C1[%d] = %0.5f, C3[%d] = %0.5f\n", N, C1[0 + N], M, C3[0 + M]);
    printf("C1[%d] = %0.5f, C3[%d] = %0.5f\n", N + 1, C1[0 + N + 1], M + 1, C3[0 + M + 1]);
    printf("\n");
    //*/

    /*******************************************
     *     matrix_fast_matmult_sse2_4x2()      *
     *******************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matrix_fast_matmult_sse2_4x2(M, K, N, A, B, C3);
    stopWatch.stop();
    elapsedTime2 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matrix_fast_matmult_sse2_4x2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#else
    printf("[matrix_fast_matmult_sse2_4x2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*****************************************
     *     matmult_s_tiling_sse2()      *
     *****************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_tiling_sse2(M, K, N, A, B, C4);
    stopWatch.stop();
    elapsedTime2 = stopWatch.getMillisec();

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_tiling_sse2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#else
    printf("[matmult_s_tiling_sse2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

    /*********************************************
     *     matmult_s_tiling_sse2_4x2()      *
     *********************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_tiling_sse2_4x2(M, K, N, A, B, C3);
    stopWatch.stop();
    elapsedTime2 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_tiling_sse2_4x2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#else
    printf("[matmult_s_tiling_sse2_4x2], elapsed time: %0.2f ms\n\n", elapsedTime2);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*****************************
     *     serial_matmult()      *
     *****************************/
#if 1
    stopWatch.start();
    serial_matmult(M, K, N, A, B, C2);
    stopWatch.stop();
    elapsedTime3 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C2, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[serial_matmult],            elapsed time: %0.2f ms\n\n", elapsedTime3);
#else
    printf("[serial_matmult],     elapsed time: %0.2f ms\n\n", elapsedTime3);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#else
    elapsedTime3 = 0.0;
#endif

    //goto MATRIX_MULT_EXIT;

#if 0
    /*******************************************
     *     matmult_s_row_MxN_K()     *
     *******************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_MxN_K(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_MxN_K], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_MxN_K], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */
#endif  /* matmult_s_row_MxN_K() */

#if 0
    /*******************************************
     *     matmult_s_row_NxM_K()     *
     *******************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_NxM_K(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_NxM_K], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_NxM_K], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */
#endif  /* matmult_s_row_NxM_K() */

    /**************************************************
     *     matmult_s_row_MxN_K_transB()     *
     **************************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_MxN_K_transB(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_MxN_K_transB], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_MxN_K_transB], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /**************************************************
     *     matmult_s_row_NxM_K_transB()     *
     **************************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_NxM_K_transB(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_NxM_K_transB], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_NxM_K_transB], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /**********************************************************
     *     matmult_s_row_tiling_MxN_K_transB()      *
     **********************************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_MxN_K_transB(M, K, N, A, B, C3);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_tiling_MxN_K_transB], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_tiling_MxN_K_transB], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*******************************************
     *     matmult_s_row_MxK_N()     *
     *******************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_MxK_N(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_MxK_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_MxK_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*******************************************
     *     matmult_s_row_KxM_N()     *
     *******************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_KxM_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_KxM_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /**************************************************
     *     matmult_s_row_tiling_MxK_N()     *
     **************************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_MxK_N(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_tiling_MxK_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_tiling_MxK_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /**************************************************
     *     matmult_s_row_tiling_KxM_N()     *
     **************************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_tiling_KxM_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_tiling_KxM_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /********************************************************
     *     matmult_s_row_tiling_k_n_m_KxM_N()     *
     ********************************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_k_n_m_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_tiling_k_n_m_KxM_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_tiling_k_n_m_KxM_N], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*******************************************************
     *     matmult_s_row_tiling_KxM_N_sse2()     *
     *******************************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_KxM_N_sse2(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_tiling_KxM_N_sse2], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_tiling_KxM_N_sse2], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /***********************************************************
     *     matmult_s_row_tiling_KxM_N_sse2_2x4()     *
     ***********************************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_tiling_KxM_N_sse2_2x4(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_tiling_KxM_N_sse2_2x4], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_tiling_KxM_N_sse2_2x4], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#ifdef DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

#if 0
    /*******************************************
     *     matmult_s_row_KxN_M()     *
     *******************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_KxN_M(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_KxN_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_KxN_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */
#endif

#if 0
    /*******************************************
     *     matmult_s_row_NxK_M()     *
     *******************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_row_NxK_M(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_row_NxK_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_row_NxK_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */
#endif

    /*******************************************
     *     matmult_s_col_KxN_M()     *
     *******************************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_col_KxN_M(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_col_KxN_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_col_KxN_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*******************************************
     *     matmult_s_col_NxK_M()     *
     *******************************************/
    matrix_init_elements(C4, M, N, MatInitZeros);

    stopWatch.start();
    matmult_s_col_NxK_M(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[matmult_s_col_NxK_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[matmult_s_col_NxK_M], elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    /*********************************
     *     matmult_see()      *
     *********************************/
    matrix_init_elements(C3, M, N, MatInitZeros);

    stopWatch.start();
    serial_matmult_sse(M, K, N, A, B, C3);
    stopWatch.stop();
    elapsedTime4 = stopWatch.getMillisec();

    diff_nums = 0;
    verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
    printf("[serial_matmult_sse],         elapsed time: %0.2f ms\n\n", elapsedTime4);
#else
    printf("[serial_matmult_sse],  elapsed time: %0.2f ms\n\n", elapsedTime4);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

#if DISPLAY_MATRIX_COMPARE
    printf("verify_ok = %d, diff_nums = %d, equal_nums = %d\n\n", verify_ok, diff_nums, M * N - diff_nums);
#endif  /* DISPLAY_MATRIX_COMPARE */

    goto MATRIX_MULT_EXIT;

#if USE_LARGE_PAGES
MATRIX_MULT_EXIT:
    if (map_address) {
        if (huge_tlb_free(map_address))
            printf("huge_tlb_free() done.\n\n");
        else
            printf("huge_tlb_free() fail.\n\n");
        map_address = NULL;
    }
#else  /* !USE_LARGE_PAGES */
MATRIX_MULT_EXIT:
    if (A)  matrix_free(A);
    if (B)  matrix_free(B);
    if (C1) matrix_free(C1);
    if (C2) matrix_free(C2);
    if (C3) matrix_free(C3);
    if (C4) matrix_free(C4);

    if (tmp1) matrix_free(tmp1);
    if (tmp2) matrix_free(tmp2);
#endif  /* USE_LARGE_PAGES */
}
