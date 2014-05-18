
#ifndef _COMMON_ASM_H_
#define _COMMON_ASM_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "fast_matmult/common.h"

#if (_WIN32 || _WIN64) && _MSC_VER

#ifndef _OS_WINDOWS_
#define _OS_WINDOWS_        1
#endif

//#include <tchar.h>
#include <intrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#endif  /* (_WIN32 || _WIN64) && _MSC_VER */

#ifndef _COMMON_H_

#ifndef _ASSEMBLER_
    #ifdef QUAD_PRECISION
        typedef struct {
            unsigned long x[2];
        } xdouble;
    #elif defined EXPRECISION
        #define xdouble long double
    #else
        #define xdouble double
    #endif

    #if defined(_OS_WINDOWS_) && (defined(__64BIT__) || defined(_WIN64))
        #ifdef _MSC_VER
            typedef __int64 BLASLONG;
            typedef unsigned __int64 BLASULONG;
        #else
            typedef long long BLASLONG;
            typedef unsigned long long BLASULONG;
        #endif
    #else
        typedef long BLASLONG;
        typedef unsigned long BLASULONG;
    #endif

    #if defined(USE_64BIT_INT) && (USE_64BIT_INT != 0)
        typedef BLASLONG blasint;
    #else
        typedef int blasint;
    #endif
#else   /* !_ASSEMBLER_ */
    #ifdef USE64BITINT
        #define INT_SHIFT       3
        #define INT_SIZE        8
    #else
        #define INT_SHIFT
        #define INT_SIZE        4
    #endif
#endif  /* _ASSEMBLER_ */

#ifdef USE_XDOUBLE
    #define FLOAT_T             xdouble
    #ifdef QUAD_PRECISION
        #define XFLOAT_T        xidouble
    #endif
    #ifdef QUAD_PRECISION
        #define FLOAT_SIZE      32
        #define BASE_SHIFT      5
        #define ZBASE_SHIFT     6
    #else
        #define FLOAT_SIZE      16
        #define BASE_SHIFT      4
        #define ZBASE_SHIFT     5
    #endif
#elif defined(USE_DOUBLE)
    #define FLOAT_T             double
    #define FLOAT_SIZE          8
    #define BASE_SHIFT          3
    #define ZBASE_SHIFT         4
#else  /* USE_FLOAT */
    #define FLOAT_T             float
    #define FLOAT_SIZE          4
    #define BASE_SHIFT          2
    #define ZBASE_SHIFT         3
#endif  /* !USE_XDOUBLE */

#ifndef XFLOAT_T
    #define XFLOAT_T            FLOAT_T
#endif

#ifndef COMPLEX
    #define COMPSIZE            1
#else
    #define COMPSIZE            2
#endif

#endif  /* _COMMON_H_ */

#ifdef _MSC_VER
#define ALIGN_4         ALIGN 4
#define ALIGN_8         ALIGN 8
#define ALIGN_16        ALIGN 16
#define ALIGN_32        ALIGN 32
#define ALIGN_64        ALIGN 64
#elif defined(__APPLE__)
#define ALIGN_4         .align 2
#define ALIGN_8         .align 3
#define ALIGN_16        .align 4
#define ALIGN_32        .align 5
#define ALIGN_64        .align 6
#elif defined(__GNUC__)
#define ALIGN_4         .align 4
#define ALIGN_8         .align 8
#define ALIGN_16        .align 16
#define ALIGN_32        .align 32
#define ALIGN_64        .align 64
#else
#define ALIGN_4
#define ALIGN_8
#define ALIGN_16
#define ALIGN_32
#define ALIGN_64
#endif

#if defined(__GNUC__)
    #define BRANCH		.byte 0x3e
    #define NOBRANCH	.byte 0x2e
    #define PADDING		.byte 0x66;
    #define HALT	    hlt
#elif defined(_MSC_VER)
    #define BRANCH		_emit 0x3e
    #define NOBRANCH	_emit 0x2e
    #define PADDING		_emit 0x66;
    #define HALT		hlt
#else
    #define BRANCH
    #define NOBRANCH
    #define PADDING     nop
    #define HALT		hlt
#endif

#ifndef xorpd
#define xorpd           xorps
#endif

#ifndef movapd
#define movapd          movaps
#endif

#ifndef __STDCALL
#define __STDCALL       __stdcall
#endif

#ifndef __CDECL
#define __CDECL         __cdecl
#endif

#ifndef __FASTCALL
#define __FASTCALL      __fastcall
#endif

/* default call: __cdecl */
#ifndef __DEFCALL
#define __DEFCALL       __CDECL
#endif

#ifdef  USE_DOUBLE
#undef  FLOAT_T
#define FLOAT_T         double
#endif  /* USE_DOUBLE */

#endif  /* _COMMON_ASM_H_  */
