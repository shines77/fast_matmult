
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _MSC_VER
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <malloc.h>
#else
// non-msvc xmm & emm define header files
// TODO: xmm and emm header files
#endif

#include <fast_matmult/fast_matmult.h>
#include <fast_matmult/common_asm.h>
#include <fast_matmult/stop_watch.h>
#include <fast_matmult/huge_tlb.h>
#include <fast_matmult/aligned_malloc.h>

#include <fast_matmult/rowmajor/matmult_s_row_X.h>
#include <fast_matmult/rowmajor/matmult_s_row_Misc.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_K.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_N.h>
#include <fast_matmult/colmajor/matmult_s_col_X.h>

#if 1
#define _aligned_malloc             iso_aligned_malloc
#define _aligned_offset_malloc      iso_aligned_offset_malloc
#define _aligned_free               iso_aligned_free
#endif

using namespace annlab;

#define MIN(a, b)       (((a) < (b)) ? (a) : (b))
#define MAX(a, b)       (((a) > (b)) ? (a) : (b))

float_t *matrix_malloc(unsigned int M, unsigned int N,
                       unsigned int alignment /* = DEFAULT_CACHELINE */)
{
    return (float_t *)::_aligned_malloc(M * N * sizeof(float_t), alignment);
}

float_t *matrix_offset_malloc(unsigned int M, unsigned int N,
                              unsigned int alignment /* = DEFAULT_CACHELINE */,
                              unsigned int offset /* = 0 */)
{
    return (float_t *)::_aligned_offset_malloc(M * N * sizeof(float_t), alignment, offset);
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

void matrix_matmult_test(int routine_mode, unsigned int M, unsigned int K, unsigned int N)
{
    unsigned int kAlignment = 128;
    stop_watch stopWatch;
    double elapsedTime, elapsedTime1, elapsedTime2;
    int err_nums = 0, diff_nums = 0;
    bool verify_ok = false;
    unsigned int test_func_mask = TEST_FUNC_MASK_NONE;

    char *verify_result[2] = { "error", "ok"   };
    char *verify_bool[2]   = { "false", "true" };

    float_t *A, *B, *C1, *C2, *C3, *C4;
    float_t *tmp1, *tmp2;;
    float_t alpha = 1.0, beta = 0.0;

#if defined(LANG_ID) && (LANG_ID != LANG_ZH_CN)
    printf("matrix_matmult_test() start...\n\n");
#else
    printf("[矩阵乘法测试程序] 开始...\n\n");
#endif
    //printf("M = %d, N = %d, K = %d\n\n", M, N, K);

    // 包含所有纯C的测试
    if (routine_mode == TEST_FUNC_PURE_C_NO_TILING
        || routine_mode == TEST_FUNC_PURE_C_ALL
        || routine_mode == TEST_FUNC_ALL_TEST) {
        test_func_mask |= TEST_FUNC_MASK_PURE_C_NO_TILING;
    }

    // 包含所有纯C(分块)的测试
    if (routine_mode == TEST_FUNC_PURE_C_TILING
        || routine_mode == TEST_FUNC_PURE_C_ALL
        || routine_mode == TEST_FUNC_ALL_TILING
        || routine_mode == TEST_FUNC_ALL_TEST) {
        test_func_mask |= TEST_FUNC_MASK_PURE_C_TILING;
    }

    // 包含所有SSEx(分块)的测试
    if (routine_mode == TEST_FUNC_SSEX_TILING
        || routine_mode == TEST_FUNC_ALL_TILING
        || routine_mode == TEST_FUNC_ALL_TEST) {
        test_func_mask |= TEST_FUNC_MASK_SSEX_TILING;
    }

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
#elif 1
    tmp1 = tmp2 = NULL;
    kAlignment = 4096;

    // matrix A
    A  = matrix_offset_malloc(M, K, kAlignment, (kAlignment - 256));
    if (A == NULL) {
        // malloc matrix A failure
        printf("matrix_offset_malloc() matrix A failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix B
    B  = matrix_offset_malloc(K, N, kAlignment, (kAlignment - 128));
    if (B == NULL) {
        // malloc matrix B failure
        printf("matrix_offset_malloc() matrix B failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C1
    C1 = matrix_offset_malloc(M, N, kAlignment);
    if (C1 == NULL) {
        // malloc matrix C1 failure
        printf("matrix_offset_malloc() matrix C1 failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C2
    C2 = matrix_offset_malloc(M, N, kAlignment);
    if (C2 == NULL) {
        // malloc matrix C2 failure
        printf("matrix_offset_malloc() matrix C2 failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C3
    C3 = matrix_offset_malloc(M, N, kAlignment);
    if (C3 == NULL) {
        // malloc matrix C failure
        printf("matrix_offset_malloc() matrix C3 failure.\n");
        goto MATRIX_MULT_EXIT;
    }

    // matrix C4
    C4 = matrix_offset_malloc(M, N, kAlignment);
    if (C4 == NULL) {
        // malloc matrix C failure
        printf("matrix_offset_malloc() matrix C4 failure.\n");
        goto MATRIX_MULT_EXIT;
    }
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

    stopWatch.start();
    //serial_matmult(M, K, N, A, B, C1);
    matmult_s_row_tiling_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C1, N);
    stopWatch.stop();
    elapsedTime = stopWatch.getMillisec();

    verify_ok = true;
    printf("[0] 用于验证的标准算法:\n\n");
    printf("[%-38s]  time: %8.2f ms, type: rowmajor\n\n", "matmult_s_row_tiling_KxM_N", elapsedTime);

    // 所有纯C的测试
    if ((test_func_mask & TEST_FUNC_MASK_PURE_C_NO_TILING) != 0) {

        printf("[1] 不使用 tiling 分块技术的纯 C/C++ 代码:\n\n");

        /*****************************
         *     serial_matmult()      *
         *****************************/
        //matrix_init_elements(C1, M, N, MatInitZeros);

        stopWatch.start();
        serial_matmult(M, K, N, A, B, C1);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        //verify_ok = matrix_compare(C1, C2, M, N, &diff_nums);
        verify_ok = true;

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "serial_matmult", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*********************************
         *     matmult_s_row_MxK_N()     *
         *********************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_MxK_N(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_MxK_N", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*********************************
         *     matmult_s_row_KxM_N()     *
         *********************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_KxM_N", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /****************************************
         *     matmult_s_row_MxN_K_transB()     *
         ****************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_MxN_K_transB(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_MxN_K_transB", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /****************************************
         *     matmult_s_row_NxM_K_transB()     *
         ****************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_NxM_K_transB(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_NxM_K_transB", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*********************************
         *     matmult_s_col_KxN_M()     *
         *********************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_col_KxN_M(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_col_KxN_M", elapsedTime, verify_result[verify_ok]);
#if 0
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);
#endif

        /*********************************
         *     matmult_s_col_NxK_M()     *
         *********************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_col_NxK_M(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_col_NxK_M", elapsedTime, verify_result[verify_ok]);
#if 0
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);
#endif

#if 0
        /*********************************
         *     matmult_s_row_MxN_K()     *
         *********************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_MxN_K(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_MxN_K", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

#endif  /* matmult_s_row_MxN_K() */

#if 1
        /*********************************
         *     matmult_s_row_NxM_K()     *
         *********************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_NxM_K(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_NxM_K", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

#endif  /* matmult_s_row_NxM_K() */

#if 0
        /*********************************
         *     matmult_s_row_KxN_M()     *
         *********************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_KxN_M(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_KxN_M", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

#endif  /* matmult_s_row_KxN_M() */

#if 0
        /*********************************
         *     matmult_s_row_NxK_M()     *
         *********************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_NxK_M(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_NxK_M", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

#endif  /* matmult_s_row_NxK_M() */
    }

    // 所有纯C(分块)的测试
    if ((test_func_mask & TEST_FUNC_MASK_PURE_C_TILING) != 0) {

        printf("[2] 使用 tiling 分块技术的纯 C/C++ 代码:\n\n");

        /****************************************
         *     matmult_s_row_tiling_MxK_N()     *
         ****************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_MxK_N(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_MxK_N", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /****************************************
         *     matmult_s_row_tiling_KxM_N()     *
         ****************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_KxM_N", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /**********************************************
         *     matmult_s_row_tiling_k_n_m_KxM_N()     *
         **********************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_k_n_m_KxM_N(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_k_n_m_KxM_N", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*****************************************
         *     matmult_s_row_tiling_NxM_K()      *
         *****************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_NxM_K(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_NxM_K", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /************************************************
         *     matmult_s_row_tiling_MxN_K_transB()      *
         ************************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_MxN_K_transB(M, K, N, A, B, C3);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_MxN_K_transB", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

#if 0
        /*****************************************************
         *     matmult_s_row_tiling_MxN_K_transB_sse2()      *
         *****************************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_MxN_K_transB_sse2(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_MxN_K_transB_sse2", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);
#endif
    }

    // 所有SSEx(分块)的测试
    if ((test_func_mask & TEST_FUNC_MASK_SSEX_TILING) != 0) {

        printf("[3] 使用 tiling 分块技术和使用 SSEx 指令优化:\n\n");

        /*********************************
         *     serial_matmult_sse()      *
         *********************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        serial_matmult_sse(M, K, N, A, B, C3);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "serial_matmult_sse", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /**********************************
         *     matrix_fast_matmult()      *
         **********************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matrix_fast_matmult(M, K, N, A, B, C3);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

#if MATRIX_AB_NEED_TRANSPOSE
        printf("[%-38s]  time: %8.2f ms\n\n", "matrix_fast_matmult: Step1", elapsedTime1);
        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matrix_fast_matmult: Step2", elapsedTime, verify_result[verify_ok]);
#else
        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matrix_fast_matmult", elapsedTime, verify_result[verify_ok]);
#endif  /* MATRIX_AB_NEED_TRANSPOSE */

        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*******************************************
         *     matrix_fast_matmult_sse2_4x2()      *
         *******************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matrix_fast_matmult_sse2_4x2(M, K, N, A, B, C3);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matrix_fast_matmult_sse2_4x2", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*****************************************
         *     matmult_s_tiling_sse2()           *
         *****************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_tiling_sse2(M, K, N, A, B, C4);
        stopWatch.stop();
        elapsedTime2 = stopWatch.getMillisec();

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_tiling_sse2", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*********************************************
         *     matmult_s_tiling_sse2_4x2()           *
         *********************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_tiling_sse2_4x2(M, K, N, A, B, C3);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_tiling_sse2_4x2", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /*************************************************
         *     matmult_s_row_tiling_KxM_N_sse2_2x4()     *
         *************************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_KxM_N_sse2_2x4(M, N, K, 1.0, A, K, B, N, 0.0, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C4, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_KxM_N_sse2_2x4", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);

        /********************************************************
         *     matmult_s_row_tiling_KxM_N_sse2_2x4_packed()     *
         ********************************************************/
        matrix_init_elements(C3, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_KxM_N_sse2_2x4_packed(M, N, K, 1.0, A, K, B, N, 0.0, C3, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_KxM_N_sse2_2x4_pk", elapsedTime, verify_result[verify_ok]);
#if 0
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);
#endif

#if 0
        /**********************************************
         *     matmult_s_row_tiling_N_sse2_4x1()      *
         **********************************************/
        matrix_init_elements(C4, M, N, MatInitZeros);

        stopWatch.start();
        matmult_s_row_tiling_N_sse2_4x1(M, N, K, &alpha, A, K, B, N, &beta, C4, N);
        stopWatch.stop();
        elapsedTime = stopWatch.getMillisec();

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "matmult_s_row_tiling_N_sse2_4x1", elapsedTime, verify_result[verify_ok]);
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);
#endif  /* matmult_s_row_tiling_N_sse2_4x1() */

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
        elapsedTime = stopWatch.getMillisec();

        // 计算完以后再转置(还原)矩阵A
        //matrix_fast_transpose_NxN((float_t *)A, K, M);
        matrix_fast_transpose_NxN((float_t *)B, N, K);

        matrix_fast_transpose_NxN((float_t *)C3, N, M);

        //C3 -= M * N;
        diff_nums = 0;
        verify_ok = matrix_compare(C1, C3, M, N, &diff_nums);

        printf("[%-38s]  time: %8.2f ms, verify: %s\n\n", "gemm_kernel_2x4_penryn", elapsedTime, verify_result[verify_ok]);
#if 0
        if (!verify_ok)
            printf("verify_ok = %5s, diff_nums = %d, equal_nums = %d\n\n", verify_bool[verify_ok], diff_nums, M * N - diff_nums);
#endif

        /*
        printf("C1[0] = %0.5f, C3[0] = %0.5f\n", C1[0], C3[0]);
        printf("C1[1] = %0.5f, C3[1] = %0.5f\n", C1[1], C3[1]);
        printf("C1[%d] = %0.5f, C3[%d] = %0.5f\n", N, C1[0 + N], M, C3[0 + M]);
        printf("C1[%d] = %0.5f, C3[%d] = %0.5f\n", N + 1, C1[0 + N + 1], M + 1, C3[0 + M + 1]);
        printf("\n");
        //*/
    }

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
