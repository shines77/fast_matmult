
#ifndef _MATMULT_BLAS_DEF_H_
#define _MATMULT_BLAS_DEF_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#ifdef _MSC_VER
#include <fast_matmult/vs/stdint.h>
#else
#include <stdint.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* 定义我们默认使用的是否是 rowmajor 存储方式, 1为是, 0为不是 */
#define CBLAS_USE_ROWMAJOR      1

/* 矩阵使用double还是float? 1为double, 其他为float */
#define CBLAS_USE_DOUBLE        1

#if defined(CBLAS_USE_DOUBLE) && (CBLAS_USE_DOUBLE != 0)
    #ifndef CBLAS_FLOAT
        #define CBLAS_FLOAT         cblas_float
    #endif
    #define FLOAT_T_EPSILON         DOUBLE_EPSINON
    #define FLOAT_T_EPSINON_TEST    DOUBLE_EPSINON_TEST
    typedef double cblas_float;
#else
    #ifndef CBLAS_FLOAT
        #define CBLAS_FLOAT         cblas_float
    #endif
    #define FLOAT_T_EPSILON         FLOAT_EPSINON
    #define FLOAT_T_EPSINON_TEST    FLOAT_EPSINON_TEST
    typedef float cblas_float;
#endif

typedef int  cblas_int;
typedef long cblas_long;

typedef struct cblas_arg_t {
    cblas_int   m;
    cblas_int   n;
    cblas_int   k;
    void       *alpha;
    void       *a;
    cblas_int   lda;
    void       *b;
    cblas_int   ldb;
    void       *beta;
    void       *c;
    cblas_int   ldc;
    int         nthreads;
} cblas_arg_t;

#define CBLAS_INDEX size_t

/* for cblas */

typedef enum Cblas_Order {
    CBlasRowMajor = 101,
    CBlasColMajor = 102
} Cblas_Order;

typedef enum Cblas_Transpose {
    CBlasNoTrans = 111,
    CBlasTrans = 112,
    CBlasConjTrans = 113,
    CBlasConjNoTrans = 114
} Cblas_Transpose;

typedef enum Cblas_Uplo {
    CBlasUpper = 121,
    CBlasLower = 122
} Cblas_Uplo;

typedef enum Cblas_Diag {
    CBlasNonUnit = 131,
    CBlasUnit = 132
} Cblas_Diag;

typedef enum Cblas_Side {
    CBlasLeft = 141,
    CBlasRight = 142
} Cblas_Side;

/* for matrix functions */

typedef enum eMatrixOrder {
    MatColMajor = 0,
    MatRowMajor,
    MatOrderMax
} eMatrixOrder;

typedef enum eMatrixTrans {
    MatNoTrans = 0,
    MatTrans,
    MatConjTrans,
    MaTransMax
} eMatrixTrans;

typedef enum eMatrixItemOrder {
    MatItemOrderAsc = 0,
    MatItemOrderDesc,
    MatItemOrderTransAsc,
    MatItemOrderTransDesc,
    MatItemOrderMax
} eMatrixItemOrder;

typedef enum eMatrixInitFcn {
    MatInitZeros = 0,
    MatInitOnes,
    MatInitRands,
    MatInitRands_Positive,
    MatInitOrder,
    MatInitSpecify,
    MatInitFcnMax
} eMatrixInitFcn;

typedef void (*cblas_func) (int mode, cblas_arg_t *args, const cblas_float *sa, const cblas_float *sb, long dummy);

typedef void (*cblas_gemm_func) (const cblas_int m, const cblas_int n, const cblas_int k,
                                 const cblas_float alpha,
                                 const cblas_float *a, const cblas_int lda,
                                 const cblas_float *b, const cblas_int ldb,
                                 const cblas_float beta,
                                 cblas_float *c, const cblas_int ldc);

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_BLAS_DEF_H_ */
