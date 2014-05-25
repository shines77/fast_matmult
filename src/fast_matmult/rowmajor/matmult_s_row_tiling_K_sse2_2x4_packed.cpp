
#include <fast_matmult/iacaMarks.h>
#include <fast_matmult/rowmajor/matmult_s_row_tiling_K_sse2_2x4_packed.h>

#define LDA     ecx
#define LDB     ebx
#define LDC     LDB

#define AA      esi
#define BB      edi
#define CC      edx

#define KK      eax

#ifndef USE_INSTRUCTION_REORDER
#define USE_INSTRUCTION_REORDER     1
#endif

#ifndef USE_PREFETCH
#define USE_PREFETCH        0
#endif

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)

#ifndef USE_PREFETCH_C
#define USE_PREFETCH_C      1
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

void matmult_s_row_tiling_NxM_K_sse2_2x4_packed(const int M, const int N, const int K,
                                                const float_t alpha,
                                                const float_t *A, const int lda,
                                                const float_t *B, const int ldb,
                                                const float_t beta,
                                                float_t *C, const int ldc)
{
    int m, n, k;
    int m_start, m_end;
    int n_start, n_end;
    int k_start, k_end;
    int m_max, n_max, k_max;
    int m_step, n_step, k_step;

    const float_t *A_, *B_;
    float_t *C_;

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)
    m_step = 8;
    k_step = 128;
    n_step = 512;

    ///*
    // 只打开 A(t0)的预读, 317.xx ms
    m_step = 128;
    k_step = 512;
    n_step = 512;
    //*/

    ///*
    // 只打开C(t0), A(t0)的预读, 318.xx ms
    m_step = 128;
    k_step = 512;
    n_step = 512;
    //*/

    /*
    // 只打开C(t0), A(t0)的预读, 318.xx ms
    m_step = 256;
    k_step = 512;
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

    /*
    // 322.xx ms
    m_step = 256;
    k_step = 128;
    n_step = 512;
    //*/

    /*
    // 313.xx ms
    m_step = 128;
    k_step = 256;
    n_step = 512;
    //*/

    ///*
    // 310.xx ms
    m_step = 128;
    k_step = 512;
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
    //TILING_OUTER_LOOP_BEGIN(n, m, k);
    TILING_OUTER_LOOP_BEGIN(m, n, k);

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

#elif 0
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
                mov         BB, B_
                mov         CC, C_

                mov         LDA, lda
                mov         LDB, ldb

                lea         LDA, [LDA * FLOAT_SIZE]
                lea         LDB, [LDB * FLOAT_SIZE]

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

                xorpd       xmm4, xmm4
                PREFETCH_C  byte ptr [CC + 0    + 1 * FLOAT_SIZE]
                xorpd       xmm5, xmm5
                PREFETCH_C  byte ptr [CC + LDC  + 3 * FLOAT_SIZE]
                xorpd       xmm6, xmm6
                xorpd       xmm7, xmm7
#else
                xorpd       xmm4, xmm4
                xorpd       xmm5, xmm5
                xorpd       xmm6, xmm6
                xorpd       xmm7, xmm7
#endif

                mov         KK, k
                // KK = K / 8
                sar         KK, 3
                je          L15

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

                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////
#if 1
                movaps      xmm0, xmmword ptr [AA + 0 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + 2 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [BB + 0 * FLOAT_SIZE]
                movapd      xmm3, xmm2
#endif

                // CC0 = BB0 * AA0
                mulpd       xmm2, xmm0
                // CC1 = BB0 * AA1
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 2 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                // CC2 = BB1 * AA0
                mulpd       xmm2, xmm0
                // CC3 = BB1 * AA1
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 0 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + 2 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                // CC4 = BB2 * AA2
                mulpd       xmm2, xmm0
                // CC5 = BB2 * AA3
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 6 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                // CC6 = BB3 * AA2
                mulpd       xmm2, xmm0
                // CC7 = BB3 * AA3
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                ////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 4 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + 6 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [BB +  8 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 10 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 4 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + 6 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]   ; LDB * 3
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 14 * FLOAT_SIZE]   ; LDB * 3
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

    #if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
                PREFETCH_A  byte ptr [AA + (PREFETCH_SIZE_A + 8) * FLOAT_SIZE]
    #endif

                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +  8 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + 10 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [BB + 16 * FLOAT_SIZE]   ; LDB * 4
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 18 * FLOAT_SIZE]   ; LDB * 4
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +  8 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + 10 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [BB + 20 * FLOAT_SIZE]   ; LDB * 5
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 22 * FLOAT_SIZE]   ; LDB * 5
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                ////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 12 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0x44
                movaps      xmm1, xmmword ptr [AA + 14 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0x44

                movapd      xmm2, xmmword ptr [BB + 24 * FLOAT_SIZE]   ; LDB * 6
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 26 * FLOAT_SIZE]   ; LDB * 6
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 12 * FLOAT_SIZE]
                pshufd      xmm0, xmm0, 0xee
                movaps      xmm1, xmmword ptr [AA + 14 * FLOAT_SIZE]
                pshufd      xmm1, xmm1, 0xee

                movapd      xmm2, xmmword ptr [BB + 28 * FLOAT_SIZE]   ; LDB * 7
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 30 * FLOAT_SIZE]   ; LDB * 7
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

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

                add         AA, 16 * FLOAT_SIZE
                add         BB, 32 * FLOAT_SIZE

                // KK -= 1;
                sub         KK, 1
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
                mov         BB, B_
                mov         CC, C_

                mov         LDA, lda
                mov         LDB, ldb

                lea         LDA, [LDA * FLOAT_SIZE]
                lea         LDB, [LDB * FLOAT_SIZE]

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

                xorpd       xmm4, xmm4
                PREFETCH_C  byte ptr [CC + 0    + 1 * FLOAT_SIZE]
                xorpd       xmm5, xmm5
                PREFETCH_C  byte ptr [CC + LDC  + 3 * FLOAT_SIZE]
                xorpd       xmm6, xmm6
                xorpd       xmm7, xmm7
#else
                xorpd       xmm4, xmm4
                xorpd       xmm5, xmm5
                xorpd       xmm6, xmm6
                xorpd       xmm7, xmm7
#endif

                mov         KK, k
                // KK = K / 8
                sar         KK, 3
                je          L15

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

                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////
#if 1
                movaps      xmm0, xmmword ptr [AA + 0 * FLOAT_SIZE]
                punpcklqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA + 2 * FLOAT_SIZE]
                punpcklqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 0 * FLOAT_SIZE]
                movapd      xmm3, xmm2
#endif
                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 2 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 0 * FLOAT_SIZE]
                punpckhqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA + 2 * FLOAT_SIZE]
                punpckhqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 6 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                ////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +  4 * FLOAT_SIZE]
                punpcklqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA +  6 * FLOAT_SIZE]
                punpcklqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB +  8 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 10 * FLOAT_SIZE]
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +  4 * FLOAT_SIZE]
                punpckhqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA +  6 * FLOAT_SIZE]
                punpckhqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]   ; LDB * 3
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 14 * FLOAT_SIZE]   ; LDB * 3
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

    #if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
                PREFETCH_A  byte ptr [AA + (PREFETCH_SIZE_A + 8) * FLOAT_SIZE]
    #endif

                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +  8 * FLOAT_SIZE]
                punpcklqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA + 10 * FLOAT_SIZE]
                punpcklqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 16 * FLOAT_SIZE]   ; LDB * 4
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 18 * FLOAT_SIZE]   ; LDB * 4
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA +  8 * FLOAT_SIZE]
                punpckhqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA + 10 * FLOAT_SIZE]
                punpckhqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 20 * FLOAT_SIZE]   ; LDB * 5
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 22 * FLOAT_SIZE]   ; LDB * 5
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                ////////////////////////////////////////////////////

                movaps      xmm1, xmmword ptr [AA + 12 * FLOAT_SIZE]
                punpcklqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA + 14 * FLOAT_SIZE]
                punpcklqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 24 * FLOAT_SIZE]   ; LDB * 6
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 26 * FLOAT_SIZE]   ; LDB * 6
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3

                //////////////////////////////////

                movaps      xmm0, xmmword ptr [AA + 12 * FLOAT_SIZE]
                punpckhqdq  xmm0, xmm0
                movaps      xmm1, xmmword ptr [AA + 14 * FLOAT_SIZE]
                punpckhqdq  xmm1, xmm1

                movapd      xmm2, xmmword ptr [BB + 28 * FLOAT_SIZE]   ; LDB * 7
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm4, xmm2
                addpd       xmm5, xmm3

                movapd      xmm2, xmmword ptr [BB + 30 * FLOAT_SIZE]   ; LDB * 7
                movapd      xmm3, xmm2

                mulpd       xmm2, xmm0
                mulpd       xmm3, xmm1

                addpd       xmm6, xmm2
                addpd       xmm7, xmm3
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

                add         AA, 16 * FLOAT_SIZE
                add         BB, 32 * FLOAT_SIZE

                // KK -= 1;
                sub         KK, 1
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
    //TILING_OUTER_LOOP_END(n, m, k);
    TILING_OUTER_LOOP_END(m, n, k);
}
