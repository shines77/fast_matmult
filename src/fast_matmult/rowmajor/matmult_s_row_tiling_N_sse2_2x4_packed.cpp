
#include <fast_matmult/rowmajor/matmult_s_row_tiling_N_sse2_2x4_packed.h>

#define LDA     ecx
#define LDB     ebx
#define LDC     LDB

#define AA      edx
#define BB      esi
#define CC      edi

#define NN      eax

#ifndef USE_INSTRUCTION_REORDER
#define USE_INSTRUCTION_REORDER     1
#endif

#ifndef USE_PREFETCH
#define USE_PREFETCH        0
#endif

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)

#ifndef USE_PREFETCH_C
#define USE_PREFETCH_C      0
#endif

#ifndef USE_PREFETCH_A
#define USE_PREFETCH_A      1
#endif

#ifndef USE_PREFETCH_B
#define USE_PREFETCH_B      0
#endif

#else  /* !USE_PREFETCH */

#undef  USE_PREFETCH_C
#define USE_PREFETCH_C      0
#undef  USE_PREFETCH_A
#define USE_PREFETCH_A      0
#undef  USE_PREFETCH_B
#define USE_PREFETCH_B      0

#endif  /* USE_PREFETCH */

/* prefetcht0, prefetcht1, prefetchnta */

#ifndef PREFETCH_C
#define PREFETCH_C          prefetcht0
#endif

#ifndef PREFETCH_A
#define PREFETCH_A          prefetcht0
#endif

#ifndef PREFETCH_B
#define PREFETCH_B          prefetchnta
#endif

#ifndef PREFETCH_SIZE_C
#define PREFETCH_SIZE_C     (8 * 13 + 0)
#endif

#ifndef PREFETCH_SIZE_A
#define PREFETCH_SIZE_A     (8 * 15 + 4)
#endif

#ifndef PREFETCH_SIZE_B
#define PREFETCH_SIZE_B     (8 * 17 + 4)
#endif

void matmult_s_row_tiling_KxM_N_sse2_2x4_packed(const int M, const int N, const int K,
                                                const cblas_float alpha,
                                                const cblas_float *A, const int lda,
                                                const cblas_float *B, const int ldb,
                                                const cblas_float beta,
                                                cblas_float *C, const int ldc)
{
#ifdef _MSC_VER
    int m, n, k;
    int m_start, m_end;
    int n_start, n_end;
    int k_start, k_end;
    int m_max, n_max, k_max;
    int m_step, n_step, k_step;

    const cblas_float *A_, *B_;
    cblas_float *C_;
    cblas_float A_m_k = (cblas_float)0.0;

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)
    // 打开C, B的预读
    m_step = 64;
    k_step = 128;
    n_step = 512;

    /*
    // 只打开A(t0)的预读, 431.xx ms
    m_step = 128;
    k_step = 1024;
    n_step = 512;
    //*/

    /*
    // 打开A(t0), B(nta)的预读, 455.xx ms
    m_step = 128;
    k_step = 1024;
    n_step = 512;
    //*/

    /*
    // 打开A(t0), B(nta)的预读, 459.xx ms
    m_step = 64;
    k_step = 128;
    n_step = 512;
    //*/
#else
    // 435.xx ms
    m_step = 8;
    k_step = 128;
    n_step = 512;

    ///*
    // 428.xx ms
    m_step = 64;
    k_step = 128;
    n_step = 512;
    //*/
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

    // 分块外层循环顺序: K, M, N
    TILING_OUTER_LOOP_BEGIN(k, m, n);
    //TILING_OUTER_LOOP_BEGIN(m, k, n);

#if 0
        // 内层循环顺序: k, m - n
        TILING_INNER_LOOP_BEGIN_EX(k, m, 1, 1);

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

        TILING_INNER_LOOP_END_EX();
#elif 0
        // 内层循环顺序: k, m - n
        TILING_INNER_LOOP_BEGIN_EX(k, m, 4, 2);
        //TILING_INNER_LOOP_BEGIN_EX(m, k, 2, 4);

        do {
            A_ = &A[m * lda + (k << 1)];

            B_ = &B[k * ldb + (n_start << 2)];
            C_ = &C[m * ldc + n_start];

            n = (n_end - n_start);

            __asm {
                push        edi
                push        esi
                push        eax
                push        ecx
                push        ebx
                push        edx

                ///////////////////////////////////////////////////////////

                mov         AA, A_
                mov         BB, B_
                mov         CC, C_

                //mov         LDA, lda
                mov         LDB, ldb

                //lea         LDA, [LDA * FLOAT_SIZE]
                lea         LDB, [LDB * FLOAT_SIZE]

                mov         NN, n

                ALIGN_16
L01:
    #if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
                PREFETCH_C  byte ptr [CC + 0    + (PREFETCH_SIZE_C + 1) * FLOAT_SIZE]
                PREFETCH_C  byte ptr [CC + LDC  + (PREFETCH_SIZE_C + 3) * FLOAT_SIZE]
    #endif

    #if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
                PREFETCH_A  byte ptr [AA + (PREFETCH_SIZE_A + 0) * FLOAT_SIZE]
    #endif

    #if defined(USE_PREFETCH_B) && (USE_PREFETCH_B != 0)
                PREFETCH_B  byte ptr [BB + (PREFETCH_SIZE_B + 0) * FLOAT_SIZE]
    #endif

                ///////////////////////////////////////////////////////////

                movapd      xmm3, xmmword ptr [AA + 0 * FLOAT_SIZE]
                xorps       xmm4, xmm4
                pshufd      xmm0, xmm3, 0x44
                xorps       xmm5, xmm5
                pshufd      xmm1, xmm3, 0xee
                xorps       xmm6, xmm6
                movapd      xmm2, xmmword ptr [BB + 0 * FLOAT_SIZE]
                xorps       xmm7, xmm7

                mulpd       xmm2, xmm0
                addpd       xmm4, xmm2

                movapd      xmm3, xmmword ptr [BB + 2 * FLOAT_SIZE]
                mulpd       xmm3, xmm0
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
                mulpd       xmm2, xmm1
                addpd       xmm4, xmm2

                movapd      xmm3, xmmword ptr [BB + 6 * FLOAT_SIZE]
                mulpd       xmm3, xmm1
                addpd       xmm5, xmm3

                ///////////////////////////////////////////

                movapd      xmm3, xmmword ptr [AA + 2 * FLOAT_SIZE]
                pshufd      xmm0, xmm3, 0x44
                pshufd      xmm1, xmm3, 0xee

                movapd      xmm2, xmmword ptr [BB + 0 * FLOAT_SIZE]
                mulpd       xmm2, xmm0
                addpd       xmm6, xmm2

                movapd      xmm3, xmmword ptr [BB + 2 * FLOAT_SIZE]
                mulpd       xmm3, xmm0
                addpd       xmm7, xmm3

                movapd      xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
                mulpd       xmm2, xmm1
                addpd       xmm6, xmm2

                movapd      xmm3, xmmword ptr [BB + 6 * FLOAT_SIZE]
                mulpd       xmm3, xmm1
                addpd       xmm7, xmm3

                ///////////////////////////////////////////////////////////

    #if defined(USE_PREFETCH_B) && (USE_PREFETCH_B != 0)
                PREFETCH_B  byte ptr [BB + (PREFETCH_SIZE_B + 8) * FLOAT_SIZE]
    #endif

                movapd      xmm3, xmmword ptr [AA + 4 * FLOAT_SIZE]
                pshufd      xmm0, xmm3, 0x44
                pshufd      xmm1, xmm3, 0xee

                movapd      xmm2, xmmword ptr [BB + 8 * FLOAT_SIZE]
                mulpd       xmm2, xmm0
                addpd       xmm4, xmm2

                movapd      xmm3, xmmword ptr [BB + 10 * FLOAT_SIZE]
                mulpd       xmm3, xmm0
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]
                mulpd       xmm2, xmm1
                addpd       xmm4, xmm2

                movapd      xmm3, xmmword ptr [BB + 14 * FLOAT_SIZE]
                mulpd       xmm3, xmm1
                addpd       xmm5, xmm3

                ///////////////////////////////////////////

                movapd      xmm3, xmmword ptr [AA + 6 * FLOAT_SIZE]
                pshufd      xmm0, xmm3, 0x44
                pshufd      xmm1, xmm3, 0xee

                movapd      xmm2, xmmword ptr [BB + 8 * FLOAT_SIZE]
                mulpd       xmm2, xmm0
                addpd       xmm6, xmm2

                movapd      xmm3, xmmword ptr [BB + 10 * FLOAT_SIZE]
                mulpd       xmm3, xmm0
                addpd       xmm7, xmm3

                movapd      xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]
                mulpd       xmm2, xmm1
                addpd       xmm6, xmm2

                movapd      xmm3, xmmword ptr [BB + 14 * FLOAT_SIZE]
                mulpd       xmm3, xmm1
                addpd       xmm7, xmm3

                ///////////////////////////////////////////////////////////

                movapd      xmm0, xmmword ptr [CC + 0   + 0 * FLOAT_SIZE]
                movapd      xmm1, xmmword ptr [CC + 0   + 2 * FLOAT_SIZE]
                movapd      xmm2, xmmword ptr [CC + LDC + 0 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [CC + LDC + 2 * FLOAT_SIZE]

                addpd       xmm4, xmm0
                addpd       xmm5, xmm1
                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                movapd      xmmword ptr [CC + 0   + 0 * FLOAT_SIZE], xmm4
                movapd      xmmword ptr [CC + 0   + 2 * FLOAT_SIZE], xmm5
                movapd      xmmword ptr [CC + LDC + 0 * FLOAT_SIZE], xmm6
                movapd      xmmword ptr [CC + LDC + 2 * FLOAT_SIZE], xmm7

                ///////////////////////////////////////////////////////////

                add         BB, 16 * FLOAT_SIZE
                add         CC,  4 * FLOAT_SIZE

                sub         NN, 4

                jg          L01

                ///////////////////////////////////////////////////////////

                pop         edx
                pop         ebx
                pop         ecx
                pop         eax
                pop         esi
                pop         edi
            }
        } while (0);

        TILING_INNER_LOOP_END_EX();
#else
        // 内层循环顺序: k, m - n
        TILING_INNER_LOOP_BEGIN_EX(k, m, 4, 2);
        //TILING_INNER_LOOP_BEGIN_EX(m, k, 2, 4);

        do {
            A_ = &A[m * lda + (k << 1)];

            B_ = &B[k * ldb + (n_start << 2)];
            C_ = &C[m * ldc + n_start];

            n = (n_end - n_start);

            __asm {
                push        edi
                push        esi
                push        eax
                push        ecx
                push        ebx
                push        edx

                ///////////////////////////////////////////////////////////

                mov         AA, A_
                mov         BB, B_
                mov         CC, C_

                mov         LDC, ldc
                lea         LDC, [LDC * FLOAT_SIZE]

                mov         NN, n

                ALIGN_16
L01:
    #if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
                PREFETCH_C  byte ptr [CC + 0    + (PREFETCH_SIZE_C + 1) * FLOAT_SIZE]
                PREFETCH_C  byte ptr [CC + LDC  + (PREFETCH_SIZE_C + 3) * FLOAT_SIZE]
    #endif

    #if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
                PREFETCH_A  byte ptr [AA + (PREFETCH_SIZE_A + 0) * FLOAT_SIZE]
    #endif

    #if defined(USE_PREFETCH_B) && (USE_PREFETCH_B != 0)
                PREFETCH_B  byte ptr [BB + (PREFETCH_SIZE_B + 0) * FLOAT_SIZE]
    #endif

                ///////////////////////////////////////////////////////////

                movapd      xmm6, xmmword ptr [AA + 0 * FLOAT_SIZE]

                movapd      xmm2, xmmword ptr [BB + 0 * FLOAT_SIZE]
                pshufd      xmm0, xmm6, 0x44
                movapd      xmm3, xmmword ptr [BB + 2 * FLOAT_SIZE]                

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm0

                xorps       xmm4, xmm4
                xorps       xmm5, xmm5

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                pshufd      xmm1, xmm6, 0xee

                movapd      xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [BB + 6 * FLOAT_SIZE]

                mulpd       xmm2, xmm1
                mulpd       xmm3, xmm1

                movapd      xmm7, xmmword ptr [AA + 2 * FLOAT_SIZE]
                xorps       xmm6, xmm6
                pshufd      xmm0, xmm7, 0x44

                addpd       xmm4, xmm2          
                addpd       xmm5, xmm3

                ///////////////////////////////////////////

                pshufd      xmm1, xmm7, 0xee
                xorps       xmm7, xmm7

                movapd      xmm2, xmmword ptr [BB + 0 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [BB + 2 * FLOAT_SIZE]

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm0

                addpd       xmm6, xmm2               
                addpd       xmm7, xmm3

                movapd      xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [BB + 6 * FLOAT_SIZE]

                mulpd       xmm2, xmm1
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2                
                addpd       xmm7, xmm3

    #if defined(USE_PREFETCH_B) && (USE_PREFETCH_B != 0)
                PREFETCH_B  byte ptr [BB + (PREFETCH_SIZE_B + 8) * FLOAT_SIZE]
    #endif

                ///////////////////////////////////////////////////////////

                movapd      xmm3, xmmword ptr [AA + 4 * FLOAT_SIZE]
                pshufd      xmm0, xmm3, 0x44

                movapd      xmm2, xmmword ptr [BB + 8 * FLOAT_SIZE]
                pshufd      xmm1, xmm3, 0xee
                movapd      xmm3, xmmword ptr [BB + 10 * FLOAT_SIZE]

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm0

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [BB + 14 * FLOAT_SIZE]

                mulpd       xmm2, xmm1
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2 
                addpd       xmm5, xmm3

                ///////////////////////////////////////////

                movapd      xmm3, xmmword ptr [AA + 6 * FLOAT_SIZE]
                pshufd      xmm0, xmm3, 0x44

                movapd      xmm2, xmmword ptr [BB + 8 * FLOAT_SIZE]
                pshufd      xmm1, xmm3, 0xee
                movapd      xmm3, xmmword ptr [BB + 10 * FLOAT_SIZE]

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm0

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                movapd      xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [BB + 14 * FLOAT_SIZE]

                mulpd       xmm2, xmm1
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2 
                addpd       xmm7, xmm3

                ///////////////////////////////////////////////////////////

                movapd      xmm0, xmmword ptr [CC + 0   + 0 * FLOAT_SIZE]
                movapd      xmm1, xmmword ptr [CC + 0   + 2 * FLOAT_SIZE]
                movapd      xmm2, xmmword ptr [CC + LDC + 0 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [CC + LDC + 2 * FLOAT_SIZE]

                addpd       xmm4, xmm0
                addpd       xmm5, xmm1
                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                movapd      xmmword ptr [CC + 0   + 0 * FLOAT_SIZE], xmm4
                movapd      xmmword ptr [CC + 0   + 2 * FLOAT_SIZE], xmm5
                movapd      xmmword ptr [CC + LDC + 0 * FLOAT_SIZE], xmm6
                movapd      xmmword ptr [CC + LDC + 2 * FLOAT_SIZE], xmm7

                ///////////////////////////////////////////////////////////

                add         BB, 16 * FLOAT_SIZE
                add         CC,  4 * FLOAT_SIZE

                sub         NN, 4

                jg          L01

                ///////////////////////////////////////////////////////////

                pop         edx
                pop         ebx
                pop         ecx
                pop         eax
                pop         esi
                pop         edi
            }
        } while (0);

        TILING_INNER_LOOP_END_EX();
#endif

    TILING_OUTER_LOOP_END(k, m, n);
    //TILING_OUTER_LOOP_END(m, k, n);

#endif  /* _MSC_VER */
}

///////////////////////////////////////////////////////////////////////////

#ifdef  USE_INSTRUCTION_REORDER
#undef  USE_INSTRUCTION_REORDER
#define USE_INSTRUCTION_REORDER     0
#endif

#ifdef  USE_PREFETCH
#undef  USE_PREFETCH
#define USE_PREFETCH    0
#endif

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)

#ifdef  USE_PREFETCH_C
#undef  USE_PREFETCH_C
#define USE_PREFETCH_C  1
#endif

#ifdef  USE_PREFETCH_A
#undef  USE_PREFETCH_A
#define USE_PREFETCH_A  0
#endif

#ifdef  USE_PREFETCH_B
#undef  USE_PREFETCH_B
#define USE_PREFETCH_B  1
#endif

#else  /* !USE_PREFETCH */

#undef  USE_PREFETCH_C
#define USE_PREFETCH_C  0
#undef  USE_PREFETCH_A
#define USE_PREFETCH_A  0
#undef  USE_PREFETCH_B
#define USE_PREFETCH_B  0

#endif  /* USE_PREFETCH */

///////////////////////////////////////////////////////////////////////////
