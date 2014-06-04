
#include <stdio.h>

#include <fast_matmult/iacaMarks.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_K_sse2_2x4.h>
#include <fast_matmult/stop_watch.h>

using namespace annlab;

#define LDA     ecx
#define LDB     ebx
#define LDC     LDB

#define AA      esi
#define BB      edi
#define CC      edx

#define B0      edi
#define B1      eax

#define LDB_X4  CC

#define KK      eax

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

void matmult_s_row_tiling_NxM_K_sse2_2x4(const int M, const int N, const int K,
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

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)
    m_step = 8;
    k_step = 128;
    n_step = 512;

    /*
    // 只打开A(t0)的预读, 453.xx ms
    m_step = 64;
    k_step = 256;
    n_step = 512;
    //*/
#else
    // 439.xx ms
    m_step = 8;
    k_step = 128;
    n_step = 512;

    /*
    // 430.xx ms
    m_step = 64;
    k_step = 128;
    n_step = 512;
    //*/

    ///*
    // 430.xx ms
    m_step = 256;
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
    //TILING_OUTER_LOOP_BEGIN(k, m, n);
    TILING_OUTER_LOOP_BEGIN(n, m, k);
    //TILING_OUTER_LOOP_BEGIN(m, n, k);

#if 0
        // 内层循环顺序: n, m - k
        TILING_INNER_LOOP_BEGIN(n, m);

        do {
            A_ = &A[m * lda + k_start];
            B_ = &B[k_start * ldb + n];
            //C_ = &C[m * ldc + n];

            C_m_n = (float_t)0.0;

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

#elif 1
        // 内层循环顺序: n, m - k
        //TILING_INNER_LOOP_BEGIN_EX(n, m, 4, 2);
        TILING_INNER_LOOP_BEGIN_EX(m, n, 2, 4);

        do {
            A_ = &A[m * lda + k_start];
            B_ = &B[k_start * ldb + n];

            C_ = &C[m * ldc + n];

            k = (k_end - k_start);

            // for Intel Architecture Code Analyzer 2.1
            //IACA_START

            __asm {
                push        edi
                push        esi
                push        eax
                push        ecx
                push        ebx
                push        edx

                ///////////////////////////////////////////////////////////

                mov         AA, A_
                //mov         BB, B_
                //mov         CC, C_

                mov         LDA, lda
                mov         LDB, ldb

                mov         B0, B_
                //mov         B1, B_

                lea         LDA, [LDA * FLOAT_SIZE]
                lea         LDB, [LDB * FLOAT_SIZE]

                mov         LDB_X4, LDB
                //add         B1, LDB
                
                sal         LDB_X4, 2

                //mov         KK, k

//                ALIGN_16
//L01:
#if 0
                movapd      xmm2, xmmword ptr [BB +       0 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                //punpcklqdq  xmm2, xmm4
                //punpckhqdq  xmm3, xmm4

                movaps      xmm0, xmmword ptr [AA +       0 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + LDA + 0 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44
#endif

#if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
                mov         CC, C_

                // eax = CC + LDC * 2; C[2, 0]
                lea     eax, [CC + LDC * 2]

                xorpd       xmm4, xmm4
                PREFETCH_C  byte ptr [CC + 0    + 1 * FLOAT_SIZE]
                xorpd       xmm5, xmm5
                PREFETCH_C  byte ptr [CC + LDC  + 3 * FLOAT_SIZE]
                xorpd       xmm6, xmm6
                PREFETCH_C  byte ptr [eax + 0   + 1 * FLOAT_SIZE]
                xorpd       xmm7, xmm7
                PREFETCH_C  byte ptr [eax + LDC + 3 * FLOAT_SIZE]
#else
                xorpd       xmm4, xmm4
                xorpd       xmm5, xmm5
                xorpd       xmm6, xmm6
                xorpd       xmm7, xmm7
#endif

                /*
                pxor        mm0, mm0
                pxor        mm1, mm1
                pxor        mm2, mm2
                pxor        mm3, mm3
                pxor        mm4, mm4
                pxor        mm5, mm5
                pxor        mm6, mm6
                pxor        mm7, mm7
                //*/

                mov     KK, k
                // KK = K / 8
                sar     KK, 3
                mov     k, KK
                je      L15

                ALIGN_16

                // IACA_START

#ifndef IACA_MARKS_OFF

                _emit   0x0F
                _emit   0x0B

                mov     ebx, 111
                _emit   0x64
                _emit   0x67
                _emit   0x90

#endif  /* IACA_MARKS_OFF */

L12:
    #if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
                PREFETCH_A  byte ptr [AA + (PREFETCH_SIZE_A + 0) * FLOAT_SIZE]
    #endif

    #if defined(USE_PREFETCH_B) && (USE_PREFETCH_B != 0)
                PREFETCH_B  byte ptr [BB + (PREFETCH_SIZE_B + 0) * FLOAT_SIZE]
                PREFETCH_B  byte ptr [BB + LDC + (PREFETCH_SIZE_B + 0) * FLOAT_SIZE]
    #endif

    #if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
                mov         CC, C_
                PREFETCH_C  byte ptr [CC + 0    + (PREFETCH_SIZE_C + 1) * FLOAT_SIZE]
                PREFETCH_C  byte ptr [CC + LDC  + (PREFETCH_SIZE_C + 3) * FLOAT_SIZE]
    #endif

                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////

                mov         B1, B0
#if 1
                movaps      xmm0, xmmword ptr [AA +       0 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + LDA + 0 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [B0 +       0 * FLOAT_SIZE]
                movapd      xmm3, xmm2
#endif

                // CC0 = BB0 * AA0
                mulpd       xmm2, xmm0
                // CC1 = BB0 * AA1
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B0 +       2 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                // CC2 = BB1 * AA0
                mulpd       xmm2, xmm0
                // CC3 = BB1 * AA1
                mulpd       xmm3, xmm1

                add         B1, LDB

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +       0 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + LDA + 0 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [B1 +       0 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                // CC4 = BB2 * AA2
                mulpd       xmm2, xmm0
                // CC5 = BB2 * AA3
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B1 +       2 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                // CC6 = BB3 * AA2
                mulpd       xmm2, xmm0
                // CC7 = BB3 * AA3
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                ////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +       2 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + LDA + 2 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [B0 + LDB * 2 + 0 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B0 + LDB * 2 + 2 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                //add         BB, LDB

                movaps      xmm0, xmmword ptr [AA +       2 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + LDA + 2 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [B1 + LDB * 2 + 0 * FLOAT_SIZE]   ; LDB * 3
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                add         B0, LDB_X4

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B1 + LDB * 2 + 2 * FLOAT_SIZE]   ; LDB * 3
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                add         B1, LDB_X4

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //sub         BB, LDB

                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +       4 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + LDA + 4 * FLOAT_SIZE]
                //lea         BB, [BB + LDB * 4]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [B0 +       0 * FLOAT_SIZE]   ; LDB * 4
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B0 +       2 * FLOAT_SIZE]   ; LDB * 4
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +       4 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + LDA + 4 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [B1 +       0 * FLOAT_SIZE]   ; LDB * 5
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B1 +       2 * FLOAT_SIZE]   ; LDB * 5
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                ////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +       6 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + LDA + 6 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [B0 + LDB * 2 + 0 * FLOAT_SIZE]   ; LDB * 6
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [B0 + LDB * 2 + 2 * FLOAT_SIZE]   ; LDB * 6
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                //add         BB, LDB

                movaps      xmm0, xmmword ptr [AA +       6 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + LDA + 6 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [B1 + LDB * 2 + 0 * FLOAT_SIZE]   ; LDB * 7
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                //lea         CC, [LDB * 2 + LDB]

                movapd      xmm2, xmmword ptr [B1 + LDB * 2 + 2 * FLOAT_SIZE]   ; LDB * 7
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //add         BB, CC

#if 0
                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////
   
                // Prepare for next loop

                movapd      xmm2, xmmword ptr [BB +       0 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                movaps      xmm0, xmmword ptr [AA +       8 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + LDA + 8 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                ///////////////////////////////////////////////////////////////
#endif

                add         AA, 8 * FLOAT_SIZE
                add         B0, LDB_X4
                //add         B1, LDB_X4
                //add         BB, 4 * FLOAT_SIZE

                // KK -= 1;
                ///*
                mov         KK, k
                sub         KK, 1
                mov         k, KK
                //*/
                //dec         k
                BRANCH
                jne         L12

                // IACA_END

#ifndef IACA_MARKS_OFF

                mov     ebx, 222
                _emit   0x64
                _emit   0x67
                _emit   0x90

                _emit   0x0F
                _emit   0x0B

#endif  /* IACA_MARKS_OFF */

                ALIGN_16
L15:
                /*
                lea         eax, [alpha]

                // From: SSE3
                movddup     xmm3, dword ptr [eax]

                mulpd       xmm4, xmm3
                mulpd       xmm5, xmm3
                mulpd       xmm6, xmm3
                mulpd       xmm7, xmm3
                //*/

                mov         CC, C_

                ///////////////////////////////////////////////////////////

                movapd      xmm0, xmmword ptr [CC + 0   + 0 * FLOAT_SIZE]
                movapd      xmm1, xmmword ptr [CC + 0   + 2 * FLOAT_SIZE]
                movapd      xmm2, xmmword ptr [CC + LDC + 0 * FLOAT_SIZE]
                movapd      xmm3, xmmword ptr [CC + LDC + 2 * FLOAT_SIZE]

                addpd       xmm4, xmm0
                addpd       xmm6, xmm1
                addpd       xmm5, xmm2
                addpd       xmm7, xmm3

                movapd      xmmword ptr [CC + 0   + 0 * FLOAT_SIZE], xmm4
                movapd      xmmword ptr [CC + 0   + 2 * FLOAT_SIZE], xmm6
                movapd      xmmword ptr [CC + LDC + 0 * FLOAT_SIZE], xmm5
                movapd      xmmword ptr [CC + LDC + 2 * FLOAT_SIZE], xmm7

                ///////////////////////////////////////////////////////////

                // CC += 4 * FLOAT_SIZE;
                //add         CC, 4 * FLOAT_SIZE

                ///////////////////////////////////////////////////////////

                pop         edx
                pop         ebx
                pop         ecx
                pop         eax
                pop         esi
                pop         edi
            }

            // for Intel Architecture Code Analyzer 2.1
            //IACA_END

        } while (0);

        TILING_INNER_LOOP_END_EX();

#endif

    //TILING_OUTER_LOOP_END(k, m, n);
    TILING_OUTER_LOOP_END(n, m, k);
    //TILING_OUTER_LOOP_END(m, n, k);

    __asm {
        //emms
    }
}

void matmult_s_row_tiling_MxN_K_transB_sse2_2x4(const int M, const int N, const int K,
                                                const cblas_float alpha,
                                                const cblas_float *A, const int lda,
                                                const cblas_float *B, const int ldb,
                                                const cblas_float beta,
                                                cblas_float *C, const int ldc)
{
    stop_watch sw;
    int m, n, k;
    int m_start, m_end;
    int n_start, n_end;
    int k_start, k_end;
    int m_step, n_step, k_step;
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
    sw.start();
    matrix_fast_transpose_NxN((cblas_float *)B, K, N);
    sw.stop();

    printf("Transpose elapsed time: %8.2f ms\n\n", sw.getMillisec());

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
                        A_ = (cblas_float *)&A[m * K + k_start];
                        B_ = (cblas_float *)&B[n * K + k_start];
                        //C_ = &C[m * N + n];
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
                        A_ = (cblas_float *)&A[m * K + k_start];
                        B_ = (cblas_float *)&B[n * K + k_start];
                        //C_ = (float_t *)&C[m * N + n];
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
    sw.start();
    matrix_fast_transpose_NxN((cblas_float *)B, N, K);
    sw.stop();

    printf("Transpose elapsed time: %8.2f ms\n\n", sw.getMillisec());
}
