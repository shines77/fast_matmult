
#include <malloc.h>
#include <memory.h>
#include <errno.h>

#ifdef _MSC_VER
#include <crtdbg.h>
#else
#define _ASSERT(expr)  ((void)0)
#define _ASSERTE(expr) ((void)0)
#endif

#include <fast_matmult/aligned_malloc.h>

#ifndef _MSC_VER
  #ifndef _RPT1
    #define _RPT1(rptno, msg, arg1)
  #endif
#endif

static unsigned char _cAlignSignFill    = 0xE9;     /* fill no-man's sign for aligned routines */
static unsigned char _cClearSignFill    = 0x00;     /* fill no-man's sign for free routines */

/*
 * The following values are non-zero, constant, odd, large, and atypical
 *      Non-zero values help find bugs assuming zero filled data.
 *      Constant values are good so that memory filling is deterministic
 *          (to help make bugs reproducable).  Of course it is bad if
 *          the constant filling of weird values masks a bug.
 *      Mathematically odd numbers are good for finding bugs assuming a cleared
 *          lower bit.
 *      Large numbers (byte values at least) are less typical, and are good
 *          at finding bad addresses.
 *      Atypical values (i.e. not too often) are good since they typically
 *          cause early detection in code.
 *      For the case of no-man's land and free blocks, if you store to any
 *          of these locations, the memory integrity checker will detect it.
 *
 *      _bAlignLandFill has been changed from 0xBD to 0xED, to ensure that
 *      4 bytes of that (0xEDEDEDED) would give an inaccessible address under 3gb.
 */

static unsigned char s_cNoMansLandFill = 0xFD;   /* fill no-man's land with this */
static unsigned char s_cAlignLandFill  = 0xED;   /* fill no-man's land for aligned routines */
static unsigned char s_cDeadLandFill   = 0xDD;   /* fill free objects with this */
static unsigned char s_cCleanLandFill  = 0xCD;   /* fill new objects with this */

/*******************************************************************************
*static int _CheckBytes() - verify byte range set to proper value
*
*Purpose:
*       verify byte range set to proper value
*
*Entry:
*       unsigned char *pb       - pointer to start of byte range
*       unsigned char bCheck    - value byte range should be set to
*       size_t nSize            - size of byte range to be checked
*
*Return:
*       TRUE - if all bytes in range equal bcheck
*       FALSE otherwise
*
*******************************************************************************/
extern "C"
bool __cdecl _iso_check_bytes(
        unsigned char *pb,
        unsigned char bCheck,
        size_t nSize
       )
{
        bool bOkay = true;
        while (nSize--) {
            if (*pb++ != bCheck) {
/* Internal error report is just noise; calling functions all report results - JWM */
/*                _RPT3(_CRT_WARN, "memory check error at 0x%p = 0x%02X, should be 0x%02X.\n", */
/*                    (BYTE *)(pb-1),*(pb-1), bCheck); */
                bOkay = false;
            }
        }
        return bOkay;
}

//
// MS1B: BitScanForward
//
// 参考:
// http://stackoverflow.com/questions/466204/rounding-off-to-nearest-power-of-2
// http://stackoverflow.com/questions/364985/algorithm-for-finding-the-smallest-power-of-two-thats-greater-or-equal-to-a-giv
//

extern "C"
size_t __cdecl iso_next_power_of_2(size_t x)
{
#if 1
    if (x == 0)
        return 0;
    // ms1b
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
#else
    size_t ms1b = 1;
    while (ms1b < x)
        ms1b <<= 1;

    return ms1b;
#endif
}

size_t __cdecl iso_adjust_alignment(size_t alignment)
{
    alignment = iso_next_power_of_2(alignment);
    _ASSERT(IS_POWER_OF_2(alignment));

#if 1
    alignment = (alignment > sizeof(uintptr_t)) ? alignment : sizeof(uintptr_t);
    _ASSERT(alignment > 0);
#endif

    return alignment;
}

size_t __cdecl iso_get_alignment(size_t alignment)
{
    return iso_adjust_alignment(alignment);
}

size_t __cdecl iso_aligned_msize(void *ptr, size_t alignment, size_t offset)
{
    size_t header_size;     /* Size of the header block */
    size_t footer_size;     /* Size of the footer block */
    size_t total_size;      /* total size of the allocated block */
    size_t user_size;       /* size of the user block*/
    uintptr_t uintptr_offset;   /* keep the alignment of the data block */
                            /* after the sizeof(void*) aligned pointer */
                            /* to the beginning of the allocated block */

    /* HEADER SIZE + FOOTER SIZE = uintptr_offset + ALIGNMENT + SIZE OF A POINTER*/
    /* HEADER SIZE + USER SIZE + FOOTER SIZE = TOTAL SIZE */
    _ASSERT(ptr != NULL);

    ALIGN_BLOCK_HEADER *pBlockHdr = NULL; /* points to the beginning of the allocated block*/
    pBlockHdr = (ALIGN_BLOCK_HEADER *)((uintptr_t)ptr & ~(sizeof(uintptr_t) - 1)) - 1;

#ifdef __linux__
    total_size = malloc_usable_size(pBlockHdr->pvAlloc);
#else
    total_size = _msize(pBlockHdr->pvAlloc);
#endif
    header_size = (uintptr_t)ptr - (uintptr_t)(pBlockHdr->pvAlloc);
    uintptr_offset = (0 - offset) & (sizeof(uintptr_t) - 1);

    /* The align cannot be smaller than the sizeof(uintptr_t) */
    alignment = (alignment > sizeof(uintptr_t) ? alignment : sizeof(uintptr_t)) -1;
    footer_size = alignment + sizeof(ALIGN_BLOCK_HEADER) + uintptr_offset - header_size;
    user_size = total_size - header_size - footer_size;
    return user_size;
}

void * __cdecl iso_aligned_malloc(size_t size, size_t alignment)
{
    size_t alignment_mask;
    size_t alloc_size;
    uintptr_t pvAlloc, pvData;
    ALIGN_BLOCK_HEADER *pBlockHdr;
#ifdef _DEBUG
    uintptr_t nFrontPaddedSize;
    uintptr_t nLastPaddedSize;
#endif

    /* 如果不是2的幂次方, 则调整到最接近的2的幂次方 */
    _ASSERT(IS_POWER_OF_2(alignment));
#if 1
    if (NOT_IS_POWER_OF_2(alignment)) {
        alignment = iso_next_power_of_2(alignment);
        _ASSERT(IS_POWER_OF_2(alignment));
    }
#endif

    alignment = (alignment > sizeof(uintptr_t)) ? alignment : sizeof(uintptr_t);
    _ASSERT(alignment > 0);

    alignment_mask = alignment - 1;

    // alloc_size align to alignment bytes (isn't must need)
    alloc_size = size + alignment_mask + sizeof(ALIGN_BLOCK_HEADER);

    pvAlloc = (uintptr_t)::malloc(alloc_size);
    if (pvAlloc != (uintptr_t)NULL) {
        // data pointer align to alignment bytes
        pvData = (uintptr_t)((pvAlloc + alignment_mask + sizeof(ALIGN_BLOCK_HEADER))
            & (~alignment_mask));

        pBlockHdr = (ALIGN_BLOCK_HEADER *)(pvData) - 1;
        _ASSERT((uintptr_t)pBlockHdr >= pvAlloc);

#if _USE_ALIGN_SIGN_
        ::memset((void *)pBlockHdr->Sign, _cAlignSignFill, ALIGN_SIGN_SIZE);
#endif
        pBlockHdr->pvAlloc = (void *)pvAlloc;

#ifdef _DEBUG
        // for debug
        nFrontPaddedSize = (uintptr_t)pvData - (uintptr_t)pvAlloc;
        nLastPaddedSize  = (uintptr_t)alloc_size - (uintptr_t)size - nFrontPaddedSize;

        _ASSERT(nFrontPaddedSize >= sizeof(ALIGN_BLOCK_HEADER));
        _ASSERT(nLastPaddedSize >= 0);
#endif
        return (void *)pvData;
    }

    return (void *)NULL;
}

void * __cdecl iso_aligned_realloc(void *ptr, size_t new_size, size_t alignment)
{
    return iso_aligned_offset_realloc(ptr, new_size, alignment, 0);
}

void * __cdecl iso_aligned_recalloc(void *ptr, size_t count, size_t new_size, size_t alignment)
{
    return iso_aligned_offset_recalloc(ptr, count, new_size, alignment, 0);
}

void * __cdecl iso_aligned_calloc(size_t count, size_t size, size_t alignment)
{
    return iso_aligned_offset_recalloc(NULL, count, size, alignment, 0);
}

void * __cdecl iso_aligned_offset_malloc(size_t size, size_t alignment, size_t offset)
{
    size_t alignment_mask;
    size_t alloc_size;
    uintptr_t pvAlloc, pvData, uintptr_offset;
    ALIGN_BLOCK_HEADER *pBlockHdr;
#ifdef _DEBUG
    uintptr_t nFrontPaddedSize;
    uintptr_t nLastPaddedSize;
#endif

    _ASSERT(offset == 0 || offset < size);
    if (offset >= size)
        return NULL;

    /* 如果不是2的幂次方, 则调整到最接近的2的幂次方 */
    _ASSERT(IS_POWER_OF_2(alignment));
#if 1
    if (NOT_IS_POWER_OF_2(alignment)) {
        alignment = iso_next_power_of_2(alignment);
        _ASSERT(IS_POWER_OF_2(alignment));
    }
#endif

    alignment = (alignment > sizeof(uintptr_t)) ? alignment : sizeof(uintptr_t);
    _ASSERT(alignment > 0);

    alignment_mask = alignment - 1;
    uintptr_offset = (0 - offset) & (sizeof(uintptr_t) - 1);

    // alloc_size align to alignment bytes (isn't must need)
    alloc_size = size + alignment_mask + sizeof(ALIGN_BLOCK_HEADER) + uintptr_offset;

    pvAlloc = (uintptr_t)::malloc(alloc_size);
    if (pvAlloc != (uintptr_t)NULL) {
        // data pointer align to alignment bytes
        pvData = (uintptr_t)((((uintptr_t)pvAlloc + alignment_mask + sizeof(ALIGN_BLOCK_HEADER)
            + uintptr_offset + offset) & (~alignment_mask)) - offset);

        pBlockHdr = (ALIGN_BLOCK_HEADER *)(pvData - uintptr_offset) - 1;
        _ASSERT((uintptr_t)pBlockHdr >= pvAlloc);

#if _USE_ALIGN_SIGN_
        ::memset((void *)pBlockHdr->Sign, _cAlignSignFill, ALIGN_SIGN_SIZE);
#endif
        pBlockHdr->pvAlloc = (void *)pvAlloc;

#ifdef _DEBUG
        // for debug
        nFrontPaddedSize = (uintptr_t)pvData - (uintptr_t)pvAlloc;
        nLastPaddedSize  = (uintptr_t)alloc_size - (uintptr_t)size - nFrontPaddedSize;

        _ASSERT(nFrontPaddedSize >= sizeof(ALIGN_BLOCK_HEADER));
        _ASSERT(nLastPaddedSize >= 0);
#endif
        return (void *)pvData;
    }

    return (void *)NULL;
}

void * __cdecl iso_aligned_offset_realloc(void *ptr, size_t new_size, size_t alignment, size_t offset)
{
    uintptr_t pvAlloc, pvData, uintptr_offset, mov_sz;
    ALIGN_BLOCK_HEADER *pBlockHdr, *s_pBlockHdr;

    if (ptr == NULL)
        return iso_aligned_offset_malloc(new_size, alignment, offset);

    if (new_size == 0) {
        iso_aligned_free(ptr);
        return NULL;
    }

    s_pBlockHdr = (ALIGN_BLOCK_HEADER *)((uintptr_t)ptr & ~(sizeof(uintptr_t) - 1)) -1;

#ifdef _DEBUG
    if (_iso_check_bytes((unsigned char *)ptr - NO_MANS_LAND_SIZE, s_cNoMansLandFill, NO_MANS_LAND_SIZE)) {
        // We don't know where (file, linenum) ptr was allocated
        _RPT1(_CRT_ERROR, "The block at 0x%p was not allocated by _aligned routines, use realloc()", ptr);
        errno = EINVAL;
        return NULL;
    }
#endif

#if _USE_ALIGN_SIGN_
    if (!_iso_check_bytes((unsigned char *)s_pBlockHdr->Sign, _cAlignSignFill, ALIGN_SIGN_SIZE)) {
        // We don't know where (file, linenum) ptr was allocated
        _RPT1(_CRT_ERROR, "Damage before 0x%p which was allocated by aligned routine\n", ptr);
    }
#endif

#ifdef __linux__
    mov_sz = malloc_usable_size(s_pBlockHdr->pvAlloc) - ((uintptr_t)ptr - (uintptr_t)s_pBlockHdr->pvAlloc);
#else
    mov_sz = _msize(s_pBlockHdr->pvAlloc) - ((uintptr_t)ptr - (uintptr_t)s_pBlockHdr->pvAlloc);
#endif

    /* validation section */
    _ASSERT(offset == 0 || offset < new_size);

    /* 如果不是2的幂次方, 则调整到最接近的2的幂次方 */
    _ASSERT(IS_POWER_OF_2(alignment));
    if (NOT_IS_POWER_OF_2(alignment)) {
        alignment = iso_next_power_of_2(alignment);
        _ASSERT(IS_POWER_OF_2(alignment));
    }

    alignment = (alignment > sizeof(uintptr_t) ? alignment : sizeof(uintptr_t)) -1;

    uintptr_offset = (0 - offset) & (sizeof(uintptr_t) - 1);

    if ((pvAlloc = (uintptr_t)::malloc(new_size + alignment + sizeof(ALIGN_BLOCK_HEADER) + uintptr_offset)) == (uintptr_t)NULL)
        return NULL;

    pvData =((pvAlloc + alignment + sizeof(ALIGN_BLOCK_HEADER) + uintptr_offset + offset) & ~alignment) - offset;
    pBlockHdr = (ALIGN_BLOCK_HEADER *)(pvData - uintptr_offset) - 1;
#if _USE_ALIGN_SIGN_
    ::memset((void *)pBlockHdr->Sign, _cAlignSignFill, ALIGN_SIGN_SIZE);
#endif
    pBlockHdr->pvAlloc = (void *)pvAlloc;

    ::memcpy((void *)pvData, ptr, mov_sz > new_size ? new_size : mov_sz);

    ::free(s_pBlockHdr->pvAlloc);

    return (void *)pvData;
}

void * __cdecl iso_aligned_offset_recalloc(void *ptr, size_t count, size_t new_size, size_t alignment, size_t offset)
{
    size_t user_size;       /* wanted size, passed to aligned realoc */
    size_t start_fill = 0;  /* location where aligned recalloc starts to fill with 0 */
                            /* filling must start from the end of the previous user block */
                            /* and must stop at the end of the new user block */
    void * pvAlloc;         /* result of aligned recalloc*/

    /* ensure that (size * count) does not overflow */
    if (count > 0) {
        _ASSERT((_HEAP_MAXREQ / count) >= new_size);
    }

    user_size = new_size * count;

    if (ptr != NULL)
        start_fill = iso_aligned_msize(ptr, alignment, offset);

    pvAlloc = iso_aligned_offset_realloc(ptr, user_size, alignment, offset);
    if (pvAlloc != NULL) {
        if (start_fill < user_size)
            ::memset((char *)pvAlloc + start_fill, 0, user_size - start_fill);
    }
    return pvAlloc;
}

void * __cdecl iso_aligned_offset_calloc(size_t count, size_t size, size_t alignment, size_t offset)
{
    return iso_aligned_offset_recalloc(NULL, count, size, alignment, offset);
}

/* this function is unuseful */
void * __cdecl _iso_free_block_header(ALIGN_BLOCK_HEADER *pBlockHdr)
{
    void *pvAlloc = NULL;

    _ASSERT(pBlockHdr != NULL);
    if (pBlockHdr != NULL) {
        pvAlloc = pBlockHdr->pvAlloc;
        _ASSERT((void *)pBlockHdr >= pvAlloc);

        if ((void *)pBlockHdr < pvAlloc) {
            // We don't know where pvData was allocated
            _RPT1(_CRT_ERROR, "Damage before 0x%p which was allocated by aligned routine\n", pvAlloc);
            pvAlloc = (void *)-2;
        }
        else {
#if _USE_ALIGN_SIGN_
            if (_iso_check_bytes(pBlockHdr->Sign, _cAlignSignFill, ALIGN_SIGN_SIZE)) {
                // Set and fill clear sign
                ::memset(pBlockHdr->Sign, _cClearSignFill, ALIGN_SIGN_SIZE);

                // Set pvAlloc's value to NULL
                pBlockHdr->pvAlloc = NULL;

                // Free memory block if need
                if (pvAlloc != NULL) {
                    ::free(pvAlloc);
                }
            }
            else {
                // We don't know where pvData was allocated
                _RPT1(_CRT_ERROR, "Damage before 0x%p which was allocated by aligned routine\n", pvAlloc);
                pvAlloc = (void *)-1;
            }
#endif
        }
    }
    return pvAlloc;
}

void __cdecl iso_aligned_free(const void *pvData)
{
    ALIGN_BLOCK_HEADER *pBlockHdr;
    void *pvAlloc;

    _ASSERT(pvData != NULL);
    if (pvData != NULL) {
        pBlockHdr = (ALIGN_BLOCK_HEADER *)((uintptr_t)pvData & ~(sizeof(uintptr_t) - 1)) - 1;
        _ASSERT(pBlockHdr < pvData);

#ifdef _DEBUG
        if (_iso_check_bytes((unsigned char *)pvData - NO_MANS_LAND_SIZE, s_cNoMansLandFill, NO_MANS_LAND_SIZE)) {
            // We don't know where (file, linenum) pvData was allocated
            _RPT1(_CRT_ERROR, "The block at 0x%p was not allocated by _aligned routines, use free()", pvData);
            return;
        }
#endif

#if _USE_ALIGN_SIGN_
        if (!_iso_check_bytes(pBlockHdr->Sign, _cAlignSignFill, ALIGN_SIGN_SIZE)) {
            // We don't know where (file, linenum) pvData was allocated
            _RPT1(_CRT_ERROR, "Damage before 0x%p which was allocated by aligned routine\n", pvData);
        }
#endif

        pvAlloc = pBlockHdr->pvAlloc;

#ifdef _DEBUG
        if ((void *)pBlockHdr < pvAlloc) {
            // We don't know where pvData was allocated
            _RPT1(_CRT_ERROR, "Damage before 0x%p which was allocated by aligned routine\n", pvData);
        }
#endif

#ifdef _ALIGNED_FREE_CLEAR_SIGN
        // Set pvAlloc's value to NULL
        pBlockHdr->pvAlloc = NULL;

        // Set and fill clear sign
  #if _USE_ALIGN_SIGN_
        ::memset(pBlockHdr->Sign, _cClearSignFill, ALIGN_SIGN_SIZE);
  #endif
#endif  /* _ALIGNED_FREE_CLEAR_SIGN */

        // Free memory block if need
        ::free(pvAlloc);
    }
}
