
#include <stdio.h>

#include <fast_matmult/common.h>
#include <fast_matmult/cblas_dgemm.h>

/**
 *
 * 我们的 cblas_dgemm 默认是使用 rowmajor 方式的
 *
 **/
void cblas_dgemm(const Cblas_Order order, const Cblas_Transpose transA,
                 const Cblas_Transpose transB,
                 const cblas_int m, const cblas_int n, const cblas_int k,
                 const float_t *alpha,
                 const float_t *a, const cblas_int lda,
                 const float_t *b, const cblas_int ldb,
                 const float_t *beta,
                 float_t *c, const cblas_int ldc,
                 cblas_func gemm_func)
{
    cblas_arg_t args;
    int transa, transb;
    cblas_int nrowa, nrowb, info;

    float_t *buffer = NULL;
    float_t *sa, *sb;

    int mode = BLAS_DOUBLE  | BLAS_REAL;

    args.alpha = (void *)alpha;
    args.beta  = (void *)beta;
    args.nthreads = 0;

    transa = -1;
    transb = -1;
    info   =  0;

/* 注意: 我们的 cblas_dgemm 默认是使用 rowmajor 方式的, 这点跟GotoBlas2是不同的 */

#if defined(CBLAS_USE_ROWMAJOR) && (CBLAS_USE_ROWMAJOR != 0)

    /* 行优先存储模式, rowmajor */

    if (order == CBlasRowMajor) {
        args.m = m;
        args.n = n;
        args.k = k;

        args.a = (void *)a;
        args.b = (void *)b;
        args.c = (void *)c;

        args.lda = lda;
        args.ldb = ldb;
        args.ldc = ldc;

        if (transA == CBlasNoTrans) { transa = 0; }
        if (transA == CBlasTrans)   { transa = 1; }

        if (transA == CBlasConjNoTrans) { transa = 2; }
        if (transA == CBlasConjTrans)   { transa = 3; }

        if (transB == CBlasNoTrans) { transb = 0; }
        if (transB == CBlasTrans)   { transb = 1; }

        if (transB == CBlasConjNoTrans) { transb = 2; }
        if (transB == CBlasConjTrans)   { transb = 3; }

        nrowa = args.m;
        if (transa & 1) { nrowa = args.k; }
        nrowb = args.k;
        if (transb & 1) { nrowb = args.n; }

        info = -1;

        if (args.ldc < args.m) { info = 13; }
        if (args.ldb < nrowb)  { info = 10; }
        if (args.lda < nrowa)  { info =  8; }
        if (args.k < 0) { info =  5; }
        if (args.n < 0) { info =  4; }
        if (args.m < 0) { info =  3; }
        if (transb < 0) { info =  2; }
        if (transa < 0) { info =  1; }
    }

    /* 列优先存储模式, colmajor */

    if (order == CBlasColMajor) {
        args.m = n;
        args.n = m;
        args.k = k;

        args.a = (void *)b;
        args.b = (void *)a;
        args.c = (void *)c;

        args.lda = ldb;
        args.ldb = lda;
        args.ldc = ldc;

        if (transB == CBlasNoTrans) { transa = 0; }
        if (transB == CBlasTrans)   { transa = 1; }

        if (transB == CBlasConjNoTrans) { transa = 0; }
        if (transB == CBlasConjTrans)   { transa = 1; }

        if (transA == CBlasNoTrans) { transb = 0; }
        if (transA == CBlasTrans)   { transb = 1; }

        if (transA == CBlasConjNoTrans) { transb = 0; }
        if (transA == CBlasConjTrans)   { transb = 1; }

        nrowa = args.m;
        if (transa & 1) { nrowa = args.k; }
        nrowb = args.k;
        if (transb & 1) { nrowb = args.n; }

        info = -1;

        if (args.ldc < args.m) { info = 13; }
        if (args.ldb < nrowb)  { info = 10; }
        if (args.lda < nrowa)  { info =  8; }
        if (args.k < 0) { info =  5; }
        if (args.n < 0) { info =  4; }
        if (args.m < 0) { info =  3; }
        if (transb < 0) { info =  2; }
        if (transa < 0) { info =  1; }
    }

#else  /* !CBLAS_USE_ROWMAJOR */

    /* 列优先存储模式, colmajor */

    if (order == CBlasColMajor) {
        args.m = m;
        args.n = n;
        args.k = k;

        args.a = (void *)a;
        args.b = (void *)b;
        args.c = (void *)c;

        args.lda = lda;
        args.ldb = ldb;
        args.ldc = ldc;

        if (transA == CBlasNoTrans) { transa = 0; }
        if (transA == CBlasTrans)   { transa = 1; }

        if (transA == CBlasConjNoTrans) { transa = 2; }
        if (transA == CBlasConjTrans)   { transa = 3; }

        if (transB == CBlasNoTrans) { transb = 0; }
        if (transB == CBlasTrans)   { transb = 1; }

        if (transB == CBlasConjNoTrans) { transb = 2; }
        if (transB == CBlasConjTrans)   { transb = 3; }

        nrowa = args.m;
        if (transa & 1) { nrowa = args.k; }
        nrowb = args.k;
        if (transb & 1) { nrowb = args.n; }

        info = -1;

        if (args.ldc < args.m) { info = 13; }
        if (args.ldb < nrowb)  { info = 10; }
        if (args.lda < nrowa)  { info =  8; }
        if (args.k < 0) { info =  5; }
        if (args.n < 0) { info =  4; }
        if (args.m < 0) { info =  3; }
        if (transb < 0) { info =  2; }
        if (transa < 0) { info =  1; }
    }

    /* 行优先存储模式, rowmajor */

    if (order == CBlasRowMajor) {
        args.m = n;
        args.n = m;
        args.k = k;

        args.a = (void *)b;
        args.b = (void *)a;
        args.c = (void *)c;

        args.lda = ldb;
        args.ldb = lda;
        args.ldc = ldc;

        if (transB == CBlasNoTrans) { transa = 0; }
        if (transB == CBlasTrans)   { transa = 1; }

        if (transB == CBlasConjNoTrans) { transa = 0; }
        if (transB == CBlasConjTrans)   { transa = 1; }

        if (transA == CBlasNoTrans) { transb = 0; }
        if (transA == CBlasTrans)   { transb = 1; }
        if (transA == CBlasConjNoTrans) { transb = 0; }
        if (transA == CBlasConjTrans)   { transb = 1; }

        nrowa = args.m;
        if (transa & 1) { nrowa = args.k; }
        nrowb = args.k;
        if (transb & 1) { nrowb = args.n; }

        info = -1;

        if (args.ldc < args.m) { info = 13; }
        if (args.ldb < nrowb)  { info = 10; }
        if (args.lda < nrowa)  { info =  8; }
        if (args.k < 0) { info =  5; }
        if (args.n < 0) { info =  4; }
        if (args.m < 0) { info =  3; }
        if (transb < 0) { info =  2; }
        if (transa < 0) { info =  1; }
    }

#endif  /* CBLAS_USE_ROWMAJOR */

    if (info >= 0) {
        // BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
        printf("cblas_dgemm() Error: info = %d\n\n", info);
        return;
    }

    if ((args.m == 0) || (args.n == 0)) { return; }

    if (gemm_func) {
        gemm_func(mode, &args, sa, sb, 0);
    }
}
