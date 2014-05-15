
#ifndef _ASSEMBLER_
#define _ASSEMBLER_
#endif

#include "fast_matmult/common_asm.h"
#include "fast_matmult/matmult_s_row_tiling_N_sse2_4x1.h"

#define STACK           16
#define ARGS            16

#if (_WIN32 || _WIN64) && _MSC_VER

#define M               DWORD PTR [esp +  4 + STACK + ARGS]
#define N               DWORD PTR [esp +  8 + STACK + ARGS]
#define K               DWORD PTR [esp + 12 + STACK + ARGS]
#define ALPHA           DWORD PTR [esp + 16 + STACK + ARGS]
#define A               DWORD PTR [esp + 20 + STACK + ARGS]
#define ARG_LDA         DWORD PTR [esp + 24 + STACK + ARGS]
#define ARG_B           DWORD PTR [esp + 28 + STACK + ARGS]
#define ARG_LDB         DWORD PTR [esp + 32 + STACK + ARGS]
#define BETA            DWORD PTR [esp + 36 + STACK + ARGS]
#define C               DWORD PTR [esp + 40 + STACK + ARGS]
#define ARG_LDC         DWORD PTR [esp + 44 + STACK + ARGS]
#define OFFSET          DWORD PTR [esp + 48 + STACK + ARGS]

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

#ifdef NANO
#define PREFETCH_SIZE   (16 * 3 + 8)
#define PREFETCHW       prefetcht0
#define PREFETCHB       prefetcht0
#endif

#ifdef NEHALEM
#define PREFETCH_SIZE   (16 * 1 - 8)
#define PREFETCHW       prefetcht0
#define PREFETCHB       prefetcht0
#endif

#ifndef PREFETCH
#define PREFETCH        prefetcht0
#endif

#ifndef PREFETCHW
#define PREFETCHW       prefetcht0
#endif

#ifndef PREFETCHB
#define PREFETCHB       prefetcht0
#endif

#ifndef PREFETCH_SIZE
#define PREFETCH_SIZE   (16 * 13 + 8)
#endif

#if (_WIN32 || _WIN64) && _MSC_VER

__declspec(naked)
void __CDECL
matmult_s_row_tiling_N_sse2_4x1(const int m, const int n, const int k,
                                const float_t *alpha,
                                const float_t *a, const int lda,
                                const float_t *b, const int ldb,
                                const float_t *beta,
                                float_t *c, const int ldc)
{
    __asm {
        sub     esp, ARGS      // # Generate Stack Frame

        push    ebp
        push    edi
        push    esi
        push    ebx

        mov     B,   ARG_B
        mov     LDC, ARG_LDC

        jmp     L999

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
matmult_s_row_tiling_N_sse2_4x1(const int M, const int N, const int K,
                                const float_t *alpha,
                                const float_t *A, const int lda,
                                const float_t *B, const int ldb,
                                const float_t *beta,
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
