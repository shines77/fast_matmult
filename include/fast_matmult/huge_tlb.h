
#ifndef _HUGE_TLB_H_
#define _HUGE_TLB_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <windows.h>

#define HUGE_TLB_SYS_PAGESIZE   4096
#define HUGE_TLB_IS_POWER_2(v)  ((v) && (!((unsigned int)(v) & ((unsigned int)(v) - 1))))

#ifdef __cplusplus
extern "C" {
#endif

extern int s_huge_tlb_inited;
extern int s_large_pagesize;

int  huge_tlb_init();
void huge_tlb_exit();

int  huge_tlb_is_inited();

int  huge_tlb_get_pagesize();
bool huge_tlb_adjust_token_privilege(TCHAR *pszPrivilege, BOOL bEnable);

void *huge_tlb_malloc(size_t size, size_t *real_alloc_size);
void *huge_tlb_malloc_ex(void *base_address, size_t size, size_t *real_alloc_size);

bool huge_tlb_free(void *p);
bool huge_tlb_free_ex(void *p, size_t size);

#ifdef __cplusplus
}
#endif

#endif  /* _HUGE_TLB_H_ */
