
#ifndef _ISO_ALIGNED_MALLOC_H
#define _ISO_ALIGNED_MALLOC_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#ifdef _MSC_VER
#include <fast_matmult/vs/stdint.h>
#else
#include <stdint.h>
#endif

#ifdef _DEBUG
#define _USE_ALIGN_SIGN_            1
#else
#define _USE_ALIGN_SIGN_            0
#endif

#define DEFAILT_CACHE_ALIGNMENT     128

// intel 32bit CPU cache line align size is 32 or 64 byte
#define MIN_MEM_ALIGNMENT           sizeof(uintptr_t)
// <= 32M Bytes (太大了内存吃不消, 浪费太多空间)
#define MAX_MEM_ALIGNMENT           0x02000000UL

/* aligned_block的标志位大小, 等于void *指针的大小 */
#define ALIGN_SIGN_SIZE             sizeof(void *)

#define NO_MANS_LAND_SIZE           4

#define UINTPTR_T_OFFSET(offset)    ((0 - offset) & (sizeof(uintptr_t) - 1))

/* for aligned_malloc() */

#define ADDR_ALIGNED_TO(ptr, alignment)  \
    (((uintptr_t)(ptr) + ((size_t)(alignment) - 1)) & ~(((size_t)(alignment) - 1)))

#define ADDR_ALIGNED_32(ptr, alignment)  \
    (((uint32_t)(ptr) + ((uint32_t)(alignment) - 1)) & ~(((uint32_t)(alignment) - 1)))

#define ADDR_ALIGNED_64(ptr, alignment)  \
    (((uint64_t)(ptr) + ((uint64_t)(alignment) - 1)) & ~(((uint64_t)(alignment) - 1)))

/* for aligned_offset_malloc() */

#define ADDR_ALIGNED_OFFSET_TO(ptr, alignment, offset)  \
    ((((uintptr_t)(ptr) + ((size_t)(alignment) - 1) + UINTPTR_T_OFFSET(offset) + (offset) \
        + (sizeof(ALIGN_BLOCK_HEADER))) & ~(((size_t)(alignment) - 1))) - (offset))

#define ADDR_ALIGNED_OFFSET_32(ptr, alignment, offset)  \
    ((((uint32_t)(ptr) + ((uint32_t)(alignment) - 1) + UINTPTR_T_OFFSET(offset) + (offset) \
        + (sizeof(ALIGN_BLOCK_HEADER))) & ~(((uint32_t)(alignment) - 1))) - (offset))

#define ADDR_ALIGNED_OFFSET_64(ptr, alignment, offset)  \
    ((((uint64_t)(ptr) + ((uint64_t)(alignment) - 1) + UINTPTR_T_OFFSET(offset) + (offset) \
        + (sizeof(ALIGN_BLOCK_HEADER))) & ~(((uint64_t)(alignment) - 1))) - (offset))

/* unsigned int v; */
/* f = (v & (v - 1)) == 0; */
/* f = v && !(v & (v - 1)); */
/* 两者是一样的, 看喜好... */

#define IS_POWER_OF_2(v)            (((size_t)(v) & ((size_t)(v) - 1)) == 0)

#define NOT_IS_POWER_OF_2(v)        ((size_t)(v) & ((size_t)(v) - 1))

#define IS_POWER_OF_2_SAFE(v)       ((v) && (!((size_t)(v) & ((size_t)(v) - 1))))

#define _IS_POWER_OF_2_SAFE(v)      ((((v) != 0)) && (((v) & ((v) - 1)) == 0))

#ifdef __cplusplus
extern "C" {
#endif

typedef struct align_block_header_t
{
    void *          pvAlloc;
#if _USE_ALIGN_SIGN_
    unsigned char   Sign[ALIGN_SIGN_SIZE];
#endif
} ALIGN_BLOCK_HEADER, * PALIGN_BLOCK_HEADER;

#ifndef iso_aligned_malloc

size_t __cdecl iso_get_alignment(size_t alignment);
size_t __cdecl iso_adjust_alignment(size_t alignment);

size_t __cdecl iso_aligned_msize(void *ptr, size_t alignment, size_t offset);

void * __cdecl iso_aligned_malloc(size_t size, size_t alignment);
void * __cdecl iso_aligned_realloc(void *ptr, size_t new_size, size_t alignment);
void * __cdecl iso_aligned_recalloc(void *ptr, size_t count, size_t new_size, size_t alignment);
void * __cdecl iso_aligned_calloc(size_t count, size_t size, size_t alignment);

void * __cdecl iso_aligned_offset_malloc(size_t size, size_t alignment, size_t offset);
void * __cdecl iso_aligned_offset_realloc(void *ptr, size_t new_size, size_t alignment, size_t offset);
void * __cdecl iso_aligned_offset_recalloc(void *ptr, size_t count, size_t new_size, size_t alignment, size_t offset);
void * __cdecl iso_aligned_offset_calloc(size_t count, size_t size, size_t alignment, size_t offset);

void   __cdecl iso_aligned_free(const void *ptr);

#endif  /* iso_aligned_malloc */

#ifdef __cplusplus
}
#endif

#endif  /* _ISO_ALIGNED_MALLOC_H */
