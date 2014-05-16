
#ifndef _COMMON_H_
#define _COMMON_H_

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

#include <fast_matmult/lang_def.h>

#endif  /* (_WIN32 || _WIN64) && _MSC_VER */

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

#endif  /* _COMMON_H_  */
