
#ifndef _FAST_MATMULT_H_
#define _FAST_MATMULT_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <stdlib.h>

#if (_WIN32 | _WIN64) && _MSC_VER
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

#include "fast_matmult/common.h"

#define USE_LARGE_PAGES         0

/* 矩阵使用double还是float? 1为double, 其他为float */
#define ANNLAB_USE_DOUBLE       1

#define FLOAT_EPSINON           (0.00001)
#define DOUBLE_EPSINON          (0.00000001)

#define FLOAT_EPSINON_TEST      (0.08f)
#define DOUBLE_EPSINON_TEST     (0.05)

#if defined(ANNLAB_USE_DOUBLE) && (ANNLAB_USE_DOUBLE != 0)
    #ifndef FLOAT_T
        #define FLOAT_T             double
    #endif
    #define FLOAT_T_EPSILON         DOUBLE_EPSINON
    #define FLOAT_T_EPSINON_TEST    DOUBLE_EPSINON_TEST
    typedef double float_t;
#else
    #ifndef FLOAT_T
        #define FLOAT_T             float
    #endif
    #define FLOAT_T_EPSILON         FLOAT_EPSINON
    #define FLOAT_T_EPSINON_TEST    FLOAT_EPSINON_TEST
    typedef float float_t;
#endif

#define DEFAULT_CACHELINE   128
#define RAND_FIXED_SEED     (2014UL)

#if _MSC_VER
#define INLINE              __inline
#define FORCE_INLINE        __forceinline
#define RESTRICT            __restrict
#else
#define INLINE              inline
#define FORCE_INLINE        inline
#define RESTRICT            restrict
#endif

#define CBLAS_ROW_MAJOR     1

#ifdef __GNUC__
    #define likely(x)       __builtin_expect(!!(x), 1)
    #define unlikely(x)     __builtin_expect(!!(x), 0)
#else
    #define likely(x)       (x)
    #define unlikely(x)     (x)
#endif

#define SetCPUAffinityMask2(m1, m2)             ((((m2) & 1) << 1) | ((m1) & 1))
#define SetCPUAffinityMask3(m1, m2, m3)         ((((m3) & 1) << 2) | (((m2) & 1) << 1) | ((m1) & 1))
#define SetCPUAffinityMask4(m1, m2, m3, m4)     ((((m4) & 1) << 3) | (((m3) & 1) << 2) | (((m2) & 1) << 1) | ((m1) & 1))

#define SetCPUAffinityMask8(m1, m2, m3, m4, m5, m6, m7, m8) \
    ((SetCPUAffinityMask4(m5, m6, m7, m8) << 4) | SetCPUAffinityMask4(m1, m2, m3, m4))

#define Lx_start(idx)       L##idx##_start
#define Lx_end(idx)         L##idx##_end
#define Lx_max(idx)         L##idx##_max

#define TILING_OUTER_LOOP_BEGIN(L1, L2, L3)                                 \
    m_max = M;                                                              \
    n_max = N;                                                              \
    k_max = K;                                                              \
    do {                                                                    \
        L1##_start = 0;                                                     \
        L1##_end   = L1##_start + L1##_step;                                \
        if (L1##_end > L1##_max)                                            \
            L1##_end = L1##_max;                                            \
        while (L1##_start < L1##_max) {                                     \
            L2##_start = 0;                                                 \
            L2##_end   = L2##_start + L2##_step;                            \
            if (L2##_end > L2##_max)                                        \
                L2##_end = L2##_max;                                        \
            while (L2##_start < L2##_max) {                                 \
                L3##_start = 0;                                             \
                L3##_end   = L3##_start + L3##_step;                        \
                if (L3##_end >L3##_max)                                     \
                    L3##_end = L3##_max;                                    \
                while (L3##_start < L3##_max) {                             \
                    do { /* do nothing */ } while (0)

#define TILING_OUTER_LOOP_END(L1, L2, L3)                                   \
                    do { /* do nothing */ } while (0);                      \
                    L3##_start += L3##_step;                                \
                    L3##_end   += L3##_step;                                \
                    if (L3##_end > L3##_max)                                \
                        L3##_end = L3##_max;                                \
                }                                                           \
                L2##_start += L2##_step;                                    \
                L2##_end   += L2##_step;                                    \
                if (L2##_end > L2##_max)                                    \
                    L2##_end = L2##_max;                                    \
            }                                                               \
            L1##_start += L1##_step;                                        \
            L1##_end   += L1##_step;                                        \
            if (L1##_end > L1##_max)                                        \
                L1##_end = L1##_max;                                        \
        }                                                                   \
    } while (0)

#define TILING_INNER_LOOP_BEGIN(L1, L2)                                     \
    do {                                                                    \
        for (L1## = L1##_start; L1## < L1##_end; ++##L1##) {                \
            for (L2## = L2##_start; L2## < L2##_end; ++##L2##) {            \
                do { /* do nothing */ } while (0)

#define TILING_INNER_LOOP_END()                                             \
            }                                                               \
        }                                                                   \
    } while (0)

#define TILING_INNER_LOOP_BEGIN_EX(L1, L2, S1, S2)                          \
    do {                                                                    \
        for (L1## = L1##_start; L1## < L1##_end; L1## += S1) {              \
            for (L2## = L2##_start; L2## < L2##_end; L2## += S2) {          \
                do { /* do nothing */ } while (0)

#define TILING_INNER_LOOP_END_EX    TILING_INNER_LOOP_END

#ifdef __cplusplus
extern "C" {
#endif

enum eMatrixOrder {
    MatColMajor = 0,
    MatRowMajor,
    MatOrderMax
};

enum eMatrixTrans {
    MatNoTrans = 0,
    MatTrans,
    MatConjTrans,
    MaTransMax
};

enum eMatrixItemOrder {
    MatItemOrderAsc = 0,
    MatItemOrderDesc,
    MatItemOrderTransAsc,
    MatItemOrderTransDesc,
    MatItemOrderMax
};

enum eMatInitFcn {
    MatInitZeros = 0,
    MatInitOnes,
    MatInitRands,
    MatInitRands_Positive,
    MatInitOrder,
    MatInitSpecify,
    MatInitFcnMax
};

enum TestFunc_Index {
    TEST_FUNC_INDEX_FIRST = 0,
    TEST_FUNC_PURE_C_NO_TILING = 1,
    TEST_FUNC_PURE_C_TILING = 2,
    TEST_FUNC_SSEX_TILING = 3,
    TEST_FUNC_PURE_C_ALL = 4,
    TEST_FUNC_ALL_TILING = 5,
    TEST_FUNC_ALL_TEST = 6,
    TEST_FUNC_INDEX_LAST
};

enum TestFunc_Mask {
    TEST_FUNC_MASK_NONE = 0,
    TEST_FUNC_MASK_PURE_C_NO_TILING = 1,
    TEST_FUNC_MASK_PURE_C_TILING = 2,
    TEST_FUNC_MASK_SSEX_TILING = 4,
    TEST_FUNC_MASK_ALL_TILING = 6,
    TEST_FUNC_MASK_ALL_TEST = 7,
    TEST_FUNC_MASK_MAX
};

void matrix_matmult_test(int routine_mode, unsigned int M, unsigned int K, unsigned int N);

void matrix_fast_transpose_NxN(float_t *A, unsigned int M, unsigned int N);

// ========================================================================

void matrix_fast_matmult(unsigned int M, unsigned int K, unsigned int N,
                         float_t *A, float_t *B, float_t *C);

void matrix_fast_matmult_sse2_4x2(unsigned int M, unsigned int K, unsigned int N,
                                  float_t *A, float_t *B, float_t *C);

// ========================================================================

void matmult_s_tiling_sse2(unsigned int M, unsigned int K, unsigned int N,
                           float_t *A, float_t *B, float_t *C);

void matmult_s_tiling_sse2_4x2(unsigned int M, unsigned int K, unsigned int N,
                               float_t *A, float_t *B, float_t *C);

// ========================================================================

void serial_matmult(unsigned int M, unsigned int K, unsigned int N,
                    float_t *A, float_t *B, float_t *C);

void serial_matmult_sse(unsigned int M, unsigned int K, unsigned int N,
                        float_t *A, float_t *B, float_t *C);

// ========================================================================

void matmult_s_row_tiling_MxN_K_transB(unsigned int M, unsigned int K, unsigned int N,
                                       float_t *A, float_t *B, float_t *C);

#ifdef __cplusplus
}
#endif

/* 因为下面这些函数用了可选参数, 所以必须以C++的方式声明,                    */
/* 如果想要在纯C环境下使用, 自己去掉可变参数就可以了(相应的要调整调用的代码) */

float_t *matrix_malloc(unsigned int M, unsigned int N,
            unsigned int alignment = DEFAULT_CACHELINE);

float_t *matrix_offset_malloc(unsigned int M, unsigned int N,
            unsigned int alignment = DEFAULT_CACHELINE, unsigned int offset = 0);

float_t *matrix_malloc_ex(unsigned int M, unsigned int N,
            unsigned int alignment = DEFAULT_CACHELINE,
            eMatInitFcn initFcn = MatInitZeros, float_t fillValue = 0.0,
            eMatrixItemOrder order = MatItemOrderAsc);

void matrix_free(float_t *A);

void matrix_init_elements(float_t *A, unsigned int M, unsigned int N,
        eMatInitFcn initFcn = MatInitZeros, float_t fillValue = 0.0,
        eMatrixItemOrder order = MatItemOrderAsc);

bool matrix_compare(const float_t *A, const float_t *B, unsigned int M, unsigned int N,
        int *diff_nums = NULL, eMatrixItemOrder order = MatItemOrderAsc);

bool matrix_transpose_verify(float_t *A, unsigned int M, unsigned int N,
        int *err_nums = NULL, eMatrixItemOrder order = MatItemOrderAsc);

#endif  /* _FAST_MATMULT_H_ */
