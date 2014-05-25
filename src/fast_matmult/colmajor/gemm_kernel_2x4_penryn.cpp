
#ifndef _ASSEMBLER_
#define _ASSEMBLER_
#endif

#include <fast_matmult/iacaMarks.h>
#include <fast_matmult/colmajor/gemm_kernel_2x4_penryn.h>

#define TRANSA          1
#undef  TRANSA

#define STACK           16
#define ARGS            16

#if (_WIN32 || _WIN64) && _MSC_VER

#define M               DWORD PTR [esp +  4 + STACK + ARGS]
#define N               DWORD PTR [esp +  8 + STACK + ARGS]
#define K               DWORD PTR [esp + 12 + STACK + ARGS]
#define ALPHA           DWORD PTR [esp + 16 + STACK + ARGS]
#define A               DWORD PTR [esp + 24 + STACK + ARGS]
#define ARG_LDA         DWORD PTR [esp + 28 + STACK + ARGS]
#define ARG_B           DWORD PTR [esp + 32 + STACK + ARGS]
#define ARG_LDB         DWORD PTR [esp + 36 + STACK + ARGS]
#define BETA            DWORD PTR [esp + 40 + STACK + ARGS]
#define C               DWORD PTR [esp + 48 + STACK + ARGS]
#define ARG_LDC         DWORD PTR [esp + 52 + STACK + ARGS]
#define OFFSET          DWORD PTR [esp + 56 + STACK + ARGS]

/* Local ARGS = 16 bytes */
#define J               DWORD PTR [esp +  0 + STACK]
#define BX              DWORD PTR [esp +  4 + STACK]
#define KK              DWORD PTR [esp +  8 + STACK]
#define KKK             DWORD PTR [esp + 12 + STACK]

#endif  /* (_WIN32 || _WIN64) && _MSC_VER */

#define AA              edx
#define BB              ecx
#define LDC             ebp
#define B               edi
#define C1              esi
#define I               ebx

#define M_STEP          4
#define K_STEP          128
#define N_STEP          512

#ifndef USE_PREFETCH
#define USE_PREFETCH    1
#endif

#if defined(USE_PREFETCH) && (USE_PREFETCH != 0)

#ifndef USE_PREFETCH_C
#define USE_PREFETCH_C  1
#endif

#ifndef USE_PREFETCH_A
#define USE_PREFETCH_A  1
#endif

#ifndef USE_PREFETCH_B
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

#ifdef NANO
#define PREFETCH_SIZE   (8 * 3 + 4)
#define PREFETCH_C      prefetcht0
#define PREFETCH_B      prefetcht0
#endif

#ifdef NEHALEM
#define PREFETCH_SIZE   (8 * 1 - 4)
#define PREFETCH_C      prefetcht0
#define PREFETCH_B      prefetcht0
#endif

/* prefetcht0, prefetcht1, prefetchnta */

#ifndef PREFETCH_C
#define PREFETCH_C      prefetcht0
#endif

#ifndef PREFETCH_A
#define PREFETCH_A      prefetcht0
#endif

#ifndef PREFETCH_B
#define PREFETCH_B      prefetcht0
#endif

#ifndef PREFETCH_SIZE
#define PREFETCH_SIZE   (8 * 13 + 4)
#endif

#if (_WIN32 || _WIN64) && _MSC_VER

__declspec(naked)
void __CDECL
gemm_kernel_2x4_penryn(const int m, const int n, const int k,
                       const float_t alpha,
                       const float_t *a, const int lda,
                       const float_t *b, const int ldb,
                       const float_t beta,
                       float_t *c, const int ldc,
                       const int offset)
{
    __asm {
        sub     esp, ARGS      // # Generate Stack Frame

        push    ebp
        push    edi
        push    esi
        push    ebx

        mov     B,   ARG_B
        mov     LDC, ARG_LDC

#ifdef TRMMKERNEL
        mov     eax, OFFSET
  #ifndef LEFT
        neg     eax
  #endif
        mov     KK, eax
#endif  /* TRMMKERNEL */

        // A += 16 * FLOAT_SIZE;
        sub     A, -16 * FLOAT_SIZE
        // B += 16 * FLOAT_SIZE;
        sub     B, -16 * FLOAT_SIZE

        // LDC = LDC * FLOAT_SIZE;
        lea     LDC, [LDC * FLOAT_SIZE]

        // J = N / 4
        mov     eax, N
        sar     eax, 2
        mov     J, eax
        jle     L30

        ALIGN_16
L01:
#if defined(TRMMKERNEL) && defined(LEFT)
        mov     eax, OFFSET
        mov     KK, eax
#endif

        // BX = ARG_B + K * FLOAT_SIZE * 4 + 16 * FLOAT_SIZE, B[4, 16]
        mov     eax, K
        sal     eax, BASE_SHIFT + 2
        lea     eax, [eax + B]
        mov     BX, eax

        mov     C1, C
        mov     AA, A

        // I = M / 2
        mov     I, M
        sar     I, 1
        jle     L20

        // IACA_START

#ifndef IACA_MARKS_OFF

        _emit   0x0F
        _emit   0x0B

        mov     ebx, 111
        _emit   0x64
        _emit   0x67
        _emit   0x90

#endif  /* IACA_MARKS_OFF */

        ALIGN_16
L11:
#if !defined(TRMMKERNEL) || \
    (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

        mov     BB, B

#else
        mov     BB, B
        mov     eax, KK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [eax * 2 + AA]
        lea     BB, [eax * 4 + BB]
#endif

        mov     eax, BX
#if defined(USE_PREFETCH_B) && (USE_PREFETCH_B != 0)
        PREFETCH_B   byte ptr [eax - 16 * FLOAT_SIZE]
#endif
        sub     eax, -8 * FLOAT_SIZE
        // BX += 8 * FLOAT_SIZE;
        mov     BX, eax

        // eax = C1 + LDC * 2; C[2, 0]
        lea     eax, [C1 + LDC * 2]

        movaps  xmm0, xmmword ptr [AA - 16 * FLOAT_SIZE]
        xorps   xmm2, xmm2
        movaps  xmm1, xmmword ptr [BB - 16 * FLOAT_SIZE]
        xorps   xmm3, xmm3

#if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
        xorps   xmm4, xmm4
        PREFETCH_C   byte ptr [C1 + 0    + 1 * FLOAT_SIZE]
        xorps   xmm5, xmm5
        PREFETCH_C   byte ptr [C1 + LDC  + 3 * FLOAT_SIZE]
        xorps   xmm6, xmm6
        PREFETCH_C   byte ptr [eax + 0   + 1 * FLOAT_SIZE]
        xorps   xmm7, xmm7
        PREFETCH_C   byte ptr [eax + LDC + 3 * FLOAT_SIZE]
#else
        xorps   xmm4, xmm4
        xorps   xmm5, xmm5
        xorps   xmm6, xmm6
        xorps   xmm7, xmm7
#endif

#ifndef TRMMKERNEL
        mov     eax, K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        mov     eax, K
        sub     eax, KK
        mov     KKK, eax
#else
        mov     eax, KK
  #ifdef LEFT
        add     eax, 2
  #else
        add     eax, 4
  #endif
        mov     KKK, eax
#endif
        // eax = K / 8
        sar     eax, 3
        je      L15

        ALIGN_16
L12:
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 0) * FLOAT_SIZE]
#endif

        //////////////////////////////////////////////

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 14 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 12 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 10 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 8 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 12 * FLOAT_SIZE]

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 6 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 4 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 10 * FLOAT_SIZE]

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 2 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 0 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 8 * FLOAT_SIZE]

        //////////////////////////////////////////////

        addpd   xmm7, xmm3
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 8) * FLOAT_SIZE]
#endif
        movaps  xmm3, xmmword ptr [BB + 2 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB + 4 * FLOAT_SIZE]
        addpd   xmm4, xmm1
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 6 * FLOAT_SIZE]

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB + 6 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB + 8 * FLOAT_SIZE]
        addpd   xmm4, xmm1
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 4 * FLOAT_SIZE]

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB + 10 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB + 12 * FLOAT_SIZE]
        addpd   xmm4, xmm1
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 2 * FLOAT_SIZE]

        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB + 14 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB + 16 * FLOAT_SIZE]
        addpd   xmm4, xmm1
        // BB += 32 * FLOAT_SIZE;
        sub     BB, -32 * FLOAT_SIZE
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 0 * FLOAT_SIZE]

        //////////////////////////////////////////////

        // AA += 16 * FLOAT_SIZE;
        sub     AA, -16 * FLOAT_SIZE

        // eax -= 1;
        sub     eax, 1
        BRANCH
        jne     L12

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
#ifndef TRMMKERNEL
        mov     eax, K
#else
        mov     eax, KKK
#endif
        and     eax, 7
        BRANCH
        je      L18

        ALIGN_16
L16:
        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 14 * FLOAT_SIZE]
        addpd   xmm6, xmm2
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 12 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        pshufd  xmm2, xmm3, 0x4e
        mulpd   xmm3, xmm0
        mulpd   xmm2, xmm0

        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]

        // AA += 2 * FLOAT_SIZE;
        add     AA, 2 * FLOAT_SIZE
        // BB += 4 * FLOAT_SIZE;
        add     BB, 4 * FLOAT_SIZE
        dec     eax
        jg      L16

        ALIGN_16
L18:
        addpd   xmm6, xmm2
        addpd   xmm7, xmm3

        // From: SSE3
        movddup xmm3, ALPHA

        movaps  xmm0, xmm4
        movsd   xmm4, xmm5
        mulpd   xmm4, xmm3
        movsd   xmm5, xmm0
        mulpd   xmm5, xmm3

        movaps  xmm0, xmm6
        movsd   xmm6, xmm7
        mulpd   xmm6, xmm3
        movsd   xmm7, xmm0
        mulpd   xmm7, xmm3

        mov     eax, C1
        or      eax, LDC
        test    eax, 15
        NOBRANCH
        jne     L18x

        /* C1 is aligned to 16 bytes address */
        lea     eax, [C1 + LDC * 2]

#ifndef TRMMKERNEL
        movaps  xmm0, xmmword ptr [C1]
        movaps  xmm1, xmmword ptr [C1 + LDC]
        movaps  xmm2, xmmword ptr [eax]
        movaps  xmm3, xmmword ptr [eax + LDC]

        addpd   xmm4, xmm0
        addpd   xmm5, xmm1
        addpd   xmm6, xmm2
        addpd   xmm7, xmm3
#endif

        movaps  xmmword ptr [C1], xmm4
        movaps  xmmword ptr [C1 + LDC], xmm5
        movaps  xmmword ptr [eax], xmm6
        movaps  xmmword ptr [eax + LDC], xmm7

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
        mov     eax, K
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [AA + eax * 2]
        lea     BB, [BB + eax * 4]
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
        add     KK, 2
#endif

        // C1 += 2 * FLOAT_SIZE;
        add     C1, 2 * FLOAT_SIZE
        // I--;
        dec     I
        jg      L11
        jmp     L20

        ALIGN_16
L18x:
        /* C1 is un-aligned to 16 bytes address */
        lea     eax, [C1 + LDC * 2]

#ifndef TRMMKERNEL
        movups  xmm0, xmmword ptr [C1]
        movups  xmm1, xmmword ptr [C1 + LDC]
        movups  xmm2, xmmword ptr [eax]
        movups  xmm3, xmmword ptr [eax + LDC]

        addpd   xmm4, xmm0
        addpd   xmm5, xmm1
        addpd   xmm6, xmm2
        addpd   xmm7, xmm3
#endif

        movups  xmmword ptr [C1], xmm4
        movups  xmmword ptr [C1 + LDC], xmm5
        movups  xmmword ptr [eax], xmm6
        movups  xmmword ptr [eax + LDC], xmm7

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
        mov     eax, K
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [AA + eax * 2]
        lea     BB, [BB + eax * 4]
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
        add     KK, 2
#endif

        // C1 += 2 * FLOAT_SIZE;
        add     C1, 2 * FLOAT_SIZE
        // I--;
        dec     I
        jg      L11

        ALIGN_16
L20:
        mov     I, M
        test    I, 1
        jle     L29

#if !defined(TRMMKERNEL) || \
    (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

        mov     BB, B

#else
        mov     BB, B
        mov     eax, KK
        lea     eax, [eax * FLOAT_SIZE]
        add     AA, eax
        lea     BB, [BB + eax * 4]
#endif

        movaps  xmm0, xmmword ptr [AA - 16 * FLOAT_SIZE]
        xorps   xmm4, xmm4
        movaps  xmm2, xmmword ptr [BB - 16 * FLOAT_SIZE]
        xorps   xmm5, xmm5
        movaps  xmm3, xmmword ptr [BB - 14 * FLOAT_SIZE]
        xorps   xmm6, xmm6
        xorps   xmm7, xmm7

#ifndef TRMMKERNEL
        mov     eax, K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        mov     eax, K
        sub     eax, KK
        mov     KKK, eax
#else
        mov     eax, KK
  #ifdef LEFT
        add     eax, 1
  #else
        add     eax, 4
  #endif
        mov     KKK, eax
#endif
        // eax = K / 8
        sar     eax, 3
        je      L25

        ALIGN_16
L22:
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 0) * FLOAT_SIZE]
#endif

        /////////////////////////////////////////////////

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 12 * FLOAT_SIZE]
        addpd   xmm5, xmm3
        movaps  xmm3, xmmword ptr [BB - 10 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm6, xmm2
        movaps  xmm2, xmmword ptr [BB - 8 * FLOAT_SIZE]
        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 6 * FLOAT_SIZE]

        /////////////////////////////////////////////////

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 4 * FLOAT_SIZE]
        addpd   xmm5, xmm3
        movaps  xmm3, xmmword ptr [BB - 2 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 12 * FLOAT_SIZE]
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm6, xmm2
        movaps  xmm2, xmmword ptr [BB - 0 * FLOAT_SIZE]
        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB - 2 * FLOAT_SIZE]

        /////////////////////////////////////////////////

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB + 4 * FLOAT_SIZE]
        addpd   xmm5, xmm3
        movaps  xmm3, xmmword ptr [BB + 6 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 10 * FLOAT_SIZE]
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm6, xmm2
        movaps  xmm2, xmmword ptr [BB +  8 * FLOAT_SIZE]
        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB + 10 * FLOAT_SIZE]

        /////////////////////////////////////////////////

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB + 12 * FLOAT_SIZE]
        addpd   xmm5, xmm3
        movaps  xmm3, xmmword ptr [BB + 14 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 8 * FLOAT_SIZE]
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm6, xmm2
        movaps  xmm2, xmmword ptr [BB + 16 * FLOAT_SIZE]
        addpd   xmm7, xmm3
        movaps  xmm3, xmmword ptr [BB + 18 * FLOAT_SIZE]

        /////////////////////////////////////////////////

        // AA +=  8 * FLOAT_SIZE;
        sub     AA,  -8 * FLOAT_SIZE
        // BB += 32 * FLOAT_SIZE;
        sub     BB, -32 * FLOAT_SIZE

        // eax -= 1;
        sub     eax, 1
        jne     L22

        ALIGN_16
L25:
#ifndef TRMMKERNEL
        mov     eax , K
#else
        mov     eax, KKK
#endif
        and     eax, 7
        BRANCH
        je      L28

        ALIGN_16
L26:
        pshufd  xmm1, xmm0, 0x44
        movsd   xmm0, xmmword ptr [AA - 15 * FLOAT_SIZE]
        mulpd   xmm2, xmm1
        mulpd   xmm3, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 12 * FLOAT_SIZE]
        addpd   xmm5, xmm3
        movaps  xmm3, xmmword ptr [BB - 10 * FLOAT_SIZE]

        // AA += 1 * FLOAT_SIZE;
        add     AA, 1 * FLOAT_SIZE
        // BB += 4 * FLOAT_SIZE;
        add     BB, 4 * FLOAT_SIZE

        dec     eax
        jg      L26

        ALIGN_16
L28:
        // From: SEE3
        movddup xmm3, ALPHA

        addpd   xmm4, xmm6
        addpd   xmm5, xmm7

        lea     eax, [C1 + LDC * 2]

#ifndef TRMMKERNEL
        movsd   xmm0, xmmword ptr [C1 + 0 * FLOAT_SIZE]
        movhpd  xmm0, xmmword ptr [C1 + LDC + 0 * FLOAT_SIZE]

        movsd   xmm1, xmmword ptr [eax + 0 * FLOAT_SIZE]
        movhpd  xmm1, xmmword ptr [eax + LDC + 0 * FLOAT_SIZE]
#endif

        mulpd   xmm4, xmm3
        mulpd   xmm5, xmm3

#ifndef TRMMKERNEL
        addpd   xmm4, xmm0
        addpd   xmm5, xmm1
#endif

        movsd   xmmword ptr [C1 + 0   + 0 * FLOAT_SIZE], xmm4
        movhpd  xmmword ptr [C1 + LDC + 0 * FLOAT_SIZE], xmm4

        movsd   xmmword ptr [eax + 0   + 0 * FLOAT_SIZE], xmm5
        movhpd  xmmword ptr [eax + LDC + 0 * FLOAT_SIZE], xmm5

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
        mov     eax, K
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        add     AA, eax
        lea     BB, [BB + eax * 4]
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
        add     KK, 1
#endif

        ALIGN_16
L29:
#if defined(TRMMKERNEL) && !defined(LEFT)
        add     KK, 4
#endif

        mov     B, BB

        lea     eax, [LDC * 4]
        add     C, eax
        dec     J
        jg      L01

        ALIGN_16
L30:
        mov     eax, N
        test    eax, 2
        jle     L50

#if defined(TRMMKERNEL) && defined(LEFT)
        mov     eax, OFFSET
        mov     KK, eax
#endif

        mov     C1, C
        mov     AA, A

        mov     I, M
        sar     I, 1
        jle     L40

        ALIGN_16
L31:
#if !defined(TRMMKERNEL) || \
    (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

        mov     BB, B

#else
        mov     BB, B
        mov     eax, KK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [AA + eax * 2]
        lea     BB, [BB + eax * 2]
#endif

        movaps  xmm0, xmmword ptr [AA - 16 * FLOAT_SIZE]
        xorps   xmm4, xmm4
        movaps  xmm1, xmmword ptr [BB - 16 * FLOAT_SIZE]
        xorps   xmm5, xmm5

#if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
        PREFETCH_C   byte ptr [C1 + 0   + 1 * FLOAT_SIZE]
        xorps   xmm6, xmm6
        PREFETCH_C   byte ptr [C1 + LDC + 1 * FLOAT_SIZE]
        xorps   xmm7, xmm7
#else
        xorps   xmm6, xmm6
        xorps   xmm7, xmm7
#endif

#ifndef TRMMKERNEL
        mov     eax, K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        mov     eax, K
        sub     eax, KK
        mov     KKK, eax
#else
        mov     eax, KK
  #ifdef LEFT
        add     eax, 2
  #else
        add     eax, 2
  #endif
        mov     KKK, eax
#endif
        // eax = K / 8
        sar     eax, 3
        je      L35

        ALIGN_16
L32:
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 0) * FLOAT_SIZE]
#endif

        /////////////////////////////////////////////////

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 14 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 12 * FLOAT_SIZE]

        addpd   xmm7, xmm1
        movaps  xmm1, xmmword ptr [BB - 12 * FLOAT_SIZE]
        addpd   xmm6, xmm2

        /////////////////////////////////////////////////

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 10 * FLOAT_SIZE]

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 10 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 8 * FLOAT_SIZE]

        addpd   xmm7, xmm1
        movaps  xmm1, xmmword ptr [BB - 8 * FLOAT_SIZE]
        addpd   xmm6, xmm2

        /////////////////////////////////////////////////

#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 8) * FLOAT_SIZE]
#endif

        /////////////////////////////////////////////////

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 6 * FLOAT_SIZE]

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 6 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 4 * FLOAT_SIZE]

        addpd   xmm7, xmm1
        movaps  xmm1, xmmword ptr [BB - 4 * FLOAT_SIZE]
        addpd   xmm6, xmm2

        /////////////////////////////////////////////////

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 2 * FLOAT_SIZE]

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 2 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 0 * FLOAT_SIZE]

        addpd   xmm7, xmm1
        movaps  xmm1, xmmword ptr [BB - 0 * FLOAT_SIZE]
        addpd   xmm6, xmm2

        /////////////////////////////////////////////////

        sub     AA, -16 * FLOAT_SIZE
        sub     BB, -16 * FLOAT_SIZE

        sub     eax, 1
        jne     L32

        ALIGN_16
L35:
#ifndef TRMMKERNEL
        mov     eax, K
#else
        mov     eax, KKK
#endif
        and     eax, 7
        BRANCH
        je      L38

        ALIGN_16
L36:
        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm1, xmm0
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]

        addpd   xmm5, xmm1
        movaps  xmm1, xmmword ptr [BB - 14 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        sub     AA, -2 * FLOAT_SIZE
        sub     BB, -2 * FLOAT_SIZE

        dec     eax
        jg      L36

        ALIGN_16
L38:
        // From: SEE3
        movddup xmm3, ALPHA

        addpd   xmm4, xmm6
        addpd   xmm5, xmm7

        movaps  xmm0, xmm4
        movsd   xmm4, xmm5
        mulpd   xmm4, xmm3
        movsd   xmm5, xmm0
        mulpd   xmm5, xmm3

#ifndef TRMMKERNEL
        movsd   xmm0, xmmword ptr [C1 + 0   + 0 * FLOAT_SIZE]
        movhpd  xmm0, xmmword ptr [C1 + 0   + 1 * FLOAT_SIZE]

        movsd   xmm1, xmmword ptr [C1 + LDC + 0 * FLOAT_SIZE]
        movhpd  xmm1, xmmword ptr [C1 + LDC + 1 * FLOAT_SIZE]

        addpd   xmm4, xmm0
        addpd   xmm5, xmm1
#endif

        movsd   xmmword ptr [C1 + 0   + 0 * FLOAT_SIZE], xmm4
        movhpd  xmmword ptr [C1 + 0   + 1 * FLOAT_SIZE], xmm4

        movsd   xmmword ptr [C1 + LDC + 0 * FLOAT_SIZE], xmm5
        movhpd  xmmword ptr [C1 + LDC + 1 * FLOAT_SIZE], xmm5

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
        mov     eax, K
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [AA + eax * 2]
        lea     BB, [BB + eax * 2]
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
        add     KK, 2
#endif

        add     C1, 2 * FLOAT_SIZE
        dec     I
        jg      L31

        ALIGN_16
L40:
        mov     I, M
        test    I, 1
        jle     L49

#if !defined(TRMMKERNEL) || \
    (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

        mov     BB, B

#else
        mov     BB, B
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        add     AA, eax
        lea     BB, [BB + eax * 2]
#endif

        movaps  xmm0, xmmword ptr [AA - 16 * FLOAT_SIZE]
        xorps   xmm4, xmm4
        movaps  xmm2, xmmword ptr [BB - 16 * FLOAT_SIZE]
        xorps   xmm5, xmm5

#ifndef TRMMKERNEL
        mov     eax, K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        mov     eax, K
        sub     eax, KK
        mov     KKK, eax
#else
        mov     eax, KK
  #ifdef LEFT
        add     eax, 1
  #else
        add     eax, 2
  #endif
        mov     KKK, eax
#endif
        // eax = K / 8
        sar     eax, 3
        je      L45

        ALIGN_16
L42:
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 0) * FLOAT_SIZE]
#endif

        /////////////////////////////////////////////

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 14 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]
        mulpd   xmm2, xmm1

        addpd   xmm5, xmm2
        movaps  xmm2, xmmword ptr [BB - 12 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 10 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 12 * FLOAT_SIZE]
        mulpd   xmm2, xmm1

        addpd   xmm5, xmm2
        movaps  xmm2, xmmword ptr [BB -  8 * FLOAT_SIZE]

        /////////////////////////////////////////////

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB -  6 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA - 10 * FLOAT_SIZE]
        mulpd   xmm2, xmm1

        addpd   xmm5, xmm2
        movaps  xmm2, xmmword ptr [BB -  4 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0x44
        mulpd   xmm2, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB -  2 * FLOAT_SIZE]

        pshufd  xmm1, xmm0, 0xee
        movaps  xmm0, xmmword ptr [AA -  8 * FLOAT_SIZE]
        mulpd   xmm2, xmm1

        addpd   xmm5, xmm2
        movaps  xmm2, xmmword ptr [BB -  0 * FLOAT_SIZE]

        /////////////////////////////////////////////

        sub     AA,  -8 * FLOAT_SIZE
        sub     BB, -16 * FLOAT_SIZE

        sub     eax, 1
        jne     L42

        ALIGN_16
L45:
#ifndef TRMMKERNEL
        mov     eax, K
#else
        mov     eax, KKK
#endif
        and     eax, 7
        BRANCH
        je      L48

        ALIGN_16
L46:
        pshufd  xmm1, xmm0, 0x44
        movsd   xmm0, xmmword ptr [AA - 15 * FLOAT_SIZE]
        mulpd   xmm2, xmm1

        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 14 * FLOAT_SIZE]

        add     AA, 1 * FLOAT_SIZE
        add     BB, 2 * FLOAT_SIZE

        dec     eax
        jg      L46

        ALIGN_16
L48:
        movddup xmm3, ALPHA

        addpd   xmm4, xmm5

#ifndef TRMMKERNEL
        movsd   xmm0, xmmword ptr [C1 + 0   + 0 * FLOAT_SIZE]
        movhpd  xmm0, xmmword ptr [C1 + LDC + 0 * FLOAT_SIZE]
#endif

        mulpd   xmm4, xmm3

#ifndef TRMMKERNEL
        addpd   xmm4, xmm0
#endif

        movsd   xmmword ptr [C1 + 0   + 0 * FLOAT_SIZE], xmm4
        movhpd  xmmword ptr [C1 + LDC + 0 * FLOAT_SIZE], xmm4

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
        mov     eax, K
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        add     AA, eax
        lea     BB, [BB + eax * 2]
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
        add     KK, 1
#endif

        ALIGN_16
L49:
#if defined(TRMMKERNEL) && !defined(LEFT)
        add     KK, 2
#endif

        mov     B, BB

        lea     eax, [LDC * 2]
        add     C, eax

        ALIGN_16
L50:
        mov     eax, N
        test    eax, 1
        jle     L999

#if defined(TRMMKERNEL) && defined(LEFT)
        mov     eax, OFFSET
        mov     KK, eax
#endif

        mov     C1, C
        mov     AA, A

        mov     I, M
        sar     I, 1
        jle     L60

        ALIGN_16
L51:
#if !defined(TRMMKERNEL) || \
    (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

        mov     BB, B
#else
        mov     BB, B
        mov     eax, KK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [AA + eax * 2]
        add     BB, eax
#endif

        movaps  xmm0, xmmword ptr [AA - 16 * FLOAT_SIZE]
        xorps   xmm4, xmm4
        movaps  xmm1, xmmword ptr [BB - 16 * FLOAT_SIZE]
        xorps   xmm5, xmm5

#if defined(USE_PREFETCH_C) && (USE_PREFETCH_C != 0)
        PREFETCH_C   byte ptr [C1 + 1 * FLOAT_SIZE]
#endif

#ifndef TRMMKERNEL
        mov     eax, K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        mov     eax, K
        sub     eax, KK
        mov     KKK, eax
#else
        mov     eax, KK
  #ifdef LEFT
        add     eax, 2
  #else
        add     eax, 1
  #endif
        mov     KKK, eax
#endif
        // eax = K / 8
        sar     eax, 3
        je      L55

        ALIGN_16
L52:
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 0) * FLOAT_SIZE]
#endif

        ////////////////////////////////////////////////

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        movaps  xmm1, xmmword ptr [BB - 14 * FLOAT_SIZE]
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 12 * FLOAT_SIZE]
        addpd   xmm5, xmm2

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 10 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        movaps  xmm1, xmmword ptr [BB - 12 * FLOAT_SIZE]
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA -  8 * FLOAT_SIZE]
        addpd   xmm5, xmm2

        /////////////////////////////////////////////////

#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 8) * FLOAT_SIZE]
#endif

        ////////////////////////////////////////////////

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA -  6 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        movaps  xmm1, xmmword ptr [BB - 10 * FLOAT_SIZE]
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA -  4 * FLOAT_SIZE]
        addpd   xmm5, xmm2

        pshufd  xmm2, xmm1, 0x4e
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA -  2 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        pshufd  xmm2, xmm1, 0x4e
        movaps  xmm1, xmmword ptr [BB -  8 * FLOAT_SIZE]
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA -  0 * FLOAT_SIZE]
        addpd   xmm5, xmm2

        ////////////////////////////////////////////////

        sub     AA, -16 * FLOAT_SIZE
        sub     BB,  -8 * FLOAT_SIZE

        sub     eax, 1
        jne     L52

        ALIGN_16
L55:
#ifndef TRMMKERNEL
        mov     eax, K
#else
        mov     eax, KKK
#endif
        and     eax, 7
        BRANCH
        je      L58

        ALIGN_16
L56:
        pshufd  xmm2, xmm1, 0x44
        movsd   xmm1, xmmword ptr [BB - 15 * FLOAT_SIZE]
        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]
        addpd   xmm4, xmm2

        add     AA, 2 * FLOAT_SIZE
        add     BB, 1 * FLOAT_SIZE

        dec     eax
        jg      L56

        ALIGN_16
L58:
        movddup xmm3, ALPHA

        addpd   xmm4, xmm5
        mulpd   xmm4, xmm3

#ifndef TRMMKERNEL
        movsd   xmm0, xmmword ptr [C1 + 0 * FLOAT_SIZE]
        movhpd  xmm0, xmmword ptr [C1 + 1 * FLOAT_SIZE]

        addpd   xmm4, xmm0
#endif

        movsd   xmmword ptr [C1 + 0 * FLOAT_SIZE], xmm4
        movhpd  xmmword ptr [C1 + 1 * FLOAT_SIZE], xmm4

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
        mov     eax, K
        sub     eax, KKK
        lea     eax, [eax * FLOAT_SIZE]
        lea     AA, [AA + eax * 2]
        add     BB, eax
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
        add     KK, 2
#endif

        add     C1, 2 * FLOAT_SIZE
        dec     I
        jg      L51

        ALIGN_16
L60:
        mov     I, M
        test    I, 1
        jle     L999

#if !defined(TRMMKERNEL) || \
    (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

        mov     BB, B
#else
        mov     BB, B
        mov     eax, KK
        lea     eax, [eax * FLOAT_SIZE]
        add     AA, eax
        add     BB, eax
#endif

        movaps  xmm0, xmmword ptr [AA - 16 * FLOAT_SIZE]
        xorps   xmm4, xmm4
        movaps  xmm2, xmmword ptr [BB - 16 * FLOAT_SIZE]
        xorps   xmm5, xmm5

#ifndef TRMMKERNEL
        mov     eax, K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        mov     eax, K
        sub     eax, KK
        mov     KKK, eax
#else
        mov     eax, KK
  #ifdef LEFT
        add     eax, 1
  #else
        add     eax, 1
  #endif
        mov     KKK, eax
#endif
        // eax = K / 8
        sar     eax, 3
        je      L65

        ALIGN_16
L62:
#if defined(USE_PREFETCH_A) && (USE_PREFETCH_A != 0)
        PREFETCH_A   byte ptr [AA + (PREFETCH_SIZE + 0) * FLOAT_SIZE]
#endif

        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 14 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 14 * FLOAT_SIZE]

        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 12 * FLOAT_SIZE]
        addpd   xmm4, xmm5
        movaps  xmm2, xmmword ptr [BB - 12 * FLOAT_SIZE]

        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 10 * FLOAT_SIZE]
        addpd   xmm4, xmm2
        movaps  xmm2, xmmword ptr [BB - 10 * FLOAT_SIZE]

        mulpd   xmm2, xmm0
        movaps  xmm0, xmmword ptr [AA - 8 * FLOAT_SIZE]
        addpd   xmm4, xmm5
        movaps  xmm2, xmmword ptr [BB - 8 * FLOAT_SIZE]

        sub     AA, -8 * FLOAT_SIZE
        sub     BB, -8 * FLOAT_SIZE

        sub     eax, 1
        jne     L62

        ALIGN_16
L65:
#ifndef TRMMKERNEL
        mov     eax, K
#else
        mov     eax, KKK
#endif
        and     eax, 7
        BRANCH
        je      L68

        ALIGN_16
L66:
        mulsd   xmm2, xmm0
        movsd   xmm0, xmmword ptr [AA - 15 * FLOAT_SIZE]
        addsd   xmm4, xmm2
        movsd   xmm2, xmmword ptr [BB - 15 * FLOAT_SIZE]

        add     AA, 1 * FLOAT_SIZE
        add     BB, 1 * FLOAT_SIZE

        dec     eax
        jg      L66

        ALIGN_16
L68:
        movddup xmm3, ALPHA

        addpd   xmm4, xmm5

        haddpd  xmm4, xmm4

#ifndef TRMMKERNEL
        movsd   xmm0, xmmword ptr [C1 + 0 * FLOAT_SIZE]
#endif

        mulsd   xmm4, xmm3

#ifndef TRMMKERNEL
        addsd   xmm4, xmm0
#endif

        movsd   xmmword ptr [C1 + 0 * FLOAT_SIZE], xmm4

        ALIGN_16
L999:
        pop     ebx
        pop     esi
        pop     edi
        pop     ebp

        add     esp, ARGS
        ret
    }
}

#else  /* !((_WIN32 || _WIN64) && _MSC_VER) */

// non-windows routine
void __CDECL
gemm_kernel_2x4_penryn(const int M, const int N, const int K,
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
    float_t A_m_k;

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

    m_max = M;
    n_max = N;
    k_max = K;

    // 分块外层循环顺序: K, M, N
    TILING_OUTER_LOOP_BEGIN(k, m, n);

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
}

#endif  /* (_WIN32 || _WIN64) && _MSC_VER */
