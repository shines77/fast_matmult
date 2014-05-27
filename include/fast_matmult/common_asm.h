
#ifndef _MATMULT_COMMON_ASM_H_
#define _MATMULT_COMMON_ASM_H_

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

#ifdef  USE_DOUBLE
#undef  FLOAT_T
#define FLOAT_T         double
#endif  /* USE_DOUBLE */

#endif  /* _MATMULT_COMMON_ASM_H_  */
