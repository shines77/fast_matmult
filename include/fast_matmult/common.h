
#ifndef _MATMULT_COMMON_H_
#define _MATMULT_COMMON_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

/* You can choose defined USE_DOUBLE, USE_XDOUBLE or USE_FLOAT */
#ifndef USE_DOUBLE
#define USE_DOUBLE          1
#endif

#define USE_64BIT_INT       0

#if (_WIN32 || _WIN64) && _MSC_VER

#ifndef _OS_WINDOWS_
#define _OS_WINDOWS_        1
#endif

#endif  /* (_WIN32 || _WIN64) && _MSC_VER */

#include <fast_matmult/lang_def.h>
#include <fast_matmult/cblas_def.h>

#define BLAS_PREC           0x0003U
#define BLAS_SINGLE         0x0000U
#define BLAS_DOUBLE         0x0001U
#define BLAS_XDOUBLE        0x0002U
#define BLAS_REAL           0x0000U
#define BLAS_COMPLEX        0x0004U

#ifndef __MMX__
#define __MMX__
#endif

#ifndef __SSE__
#define __SSE__
#endif

#ifndef __SSE2__
#define __SSE2__
#endif

#ifndef __SSE3__
#define __SSE3__
#endif

#ifdef __cplusplus
extern "C" {
#endif

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
            typedef __int64 blaslong;
            typedef unsigned __int64 blas_ulong;
        #else
            typedef long long blaslong;
            typedef unsigned long long blas_ulong;
        #endif
    #else
        typedef long blas_long;
        typedef unsigned long blas_ulong;
    #endif

    #if defined(USE_64BIT_INT) && (USE_64BIT_INT != 0)
        typedef blas_long blas_int;
    #else
        typedef int blas_int;
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

#ifdef NEEDBUNDERSCORE
#define BLASFUNC(FUNC)          FUNC##_
#else
#define BLASFUNC(FUNC)          FUNC
#endif

/* define __cdecl, __fastcall for windows */

#if defined(_MSC_VER)

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

#else  /* !_MSC_VER */

#undef __STDCALL
#undef __CDECL
#undef __FASTCALL
#undef __DEFCALL

#define __STDCALL
#define __CDECL
#define __FASTCALL
#define __DEFCALL

#endif  /* _MSC_VER */

/* define __cdecl, __fastcall for linux */

#if defined(__GNUC__)

#ifndef _GCC_STDCALL
#define _GCC_STDCALL
#endif

#ifndef _GCC_CDECL
#define _GCC_CDECL          __attribute__ ((__cdecl__))
#endif

#ifndef _GCC_FASTCALL
#define _GCC_FASTCALL       __attribute__ ((fastcall))
#endif

/* default call: __cdecl */
#ifndef _GCC_DEFCALL
#define _GCC_DEFCALL        _GCC_CDECL
#endif

#else  /* !__GNUC__ */

#undef _GCC_STDCALL
#undef _GCC_CDECL
#undef _GCC_FASTCALL
#undef _GCC_DEFCALL

#define _GCC_STDCALL
#define _GCC_CDECL
#define _GCC_FASTCALL
#define _GCC_DEFCALL

#endif  /* __GNUC__ */

/* For Keyboard virtual code */

enum MM_VT_KEYBOARDS {
    //
    MM_VT_KEY_FIRST = 0,
    MM_VT_KEY_TAB = 8,
    MM_VT_KEY_RETURN = 13,
    MM_VT_KEY_LAST = 255
};

#ifndef _USE_GETCHAR
  #ifdef __linux__
    #define _USE_GETCHAR        1
  #else
    #define _USE_GETCHAR        0
  #endif
#endif  /* _USE_GETCHAR */

#define GETCH_DEFUALT_VALUE     (-1)
#define GETCH_EXIT_PROGRAM      (0)

#ifdef __cplusplus
}
#endif

#endif  /* _MATMULT_COMMON_H_  */
