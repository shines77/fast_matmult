
/// <comment>
///
/// From: https://blogs.oracle.com/dagastine/entry/java_se_tuning_tip_large
///
///   Only Windows Server 2003/2008 (vista, windows 7, 8) supports large page memory.
/// In order to use it, the administrator must first assign additional privilege to the
/// user who will be running the application: # select Control Panel -> Administrative Tools
/// -> Local Security Policy # select Local Policies -> User Rights Assignment
/// # double click "Lock pages in memory", add users and/or groups # reboot
/// the machine As always, every application is different and true performance
/// is always defined by each individual running their own application.
///
/// </comment>

/// <comment>
///
/// Large-Page Support
///
/// From: http://msdn.microsoft.com/en-us/library/windows/desktop/aa366720%28v=vs.85%29.aspx
///
/// To use large-page support
///
/// 1. Obtain the SeLockMemoryPrivilege privilege by calling the AdjustTokenPrivileges function.
///    For more information, see Assigning Privileges to an Account and Changing Privileges in a Token.
///
/// 2. Retrieve the minimum large-page size by calling the GetLargePageMinimum function.
///
/// 3. Include the MEM_LARGE_PAGES value when calling the VirtualAlloc function.
///    The size and alignment must be a multiple of the large-page minimum.
///
/// </comment>

/// <comment>
///
/// Creating a File Mapping Using Large Pages
///
/// From: http://msdn.microsoft.com/en-us/library/windows/desktop/aa366543%28v=vs.85%29.aspx
///
/// </comment>

#include <stdio.h>
#include <windows.h>
#include <tchar.h>

#include <fast_matmult/huge_TLB.h>

typedef int (* pfn_GetLargePageMinimum)(void);

int s_huge_tlb_inited = 0;
int s_large_pagesize  = 0;

void huge_tlb_display_error(TCHAR *pszAPI, DWORD dwError)
{
    LPVOID lpvMessageBuffer = NULL;

    ::FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |
                    FORMAT_MESSAGE_FROM_SYSTEM |
                    FORMAT_MESSAGE_IGNORE_INSERTS,
                    NULL, dwError,
                    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                    (LPTSTR)&lpvMessageBuffer, 0, NULL);

    // ... now display this string
    _tprintf(_T("ERROR: API        = %s\n"), pszAPI);
    _tprintf(_T("       error code = %d\n"), dwError);
    _tprintf(_T("       message    = %s\n"), lpvMessageBuffer);

    // Free the buffer allocated by the system
    if (lpvMessageBuffer)
        ::LocalFree(lpvMessageBuffer);

    //::ExitProcess(GetLastError());
}

int huge_tlb_is_inited()
{
    return s_huge_tlb_inited;
}

int huge_tlb_init()
{
    int inited_ok = 0;
    int large_pagesize = huge_tlb_get_pagesize();
    if (large_pagesize <= 0) {
        huge_tlb_display_error(_T("huge_tlb_get_pagesize()"), 0);
        s_huge_tlb_inited = 0;
        return false;
    }

    s_large_pagesize = large_pagesize;
    if (large_pagesize < HUGE_TLB_SYS_PAGESIZE
        || !HUGE_TLB_IS_POWER_2(large_pagesize)) {
        s_huge_tlb_inited = 0;
        return false;
    }

    if (s_huge_tlb_inited == 0) {
        bool success = true;
        //success = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, TRUE);
        //success = huge_tlb_adjust_token_privilege(_T("SeLockMemoryPrivilege"), TRUE);
        if (!success)
            huge_tlb_display_error(_T("huge_tlb_adjust_token_privilege('SeLockMemoryPrivilege') failed."), 0);
        else
            inited_ok = 1;
    }
    s_huge_tlb_inited = inited_ok;
    return inited_ok;
}

void huge_tlb_exit(int echo)
{
    if (s_huge_tlb_inited != 0) {
        bool success = true;
        //success = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, FALSE);
        //success = huge_tlb_adjust_token_privilege(_T("SeLockMemoryPrivilege"), FALSE);
        if (!success) {
            huge_tlb_display_error(_T("huge_tlb_exit() failed."), 0);
        }
        else {
            if (echo)
                printf("huge_tlb_exit() suucess.\n");
        }
        s_huge_tlb_inited = 0;
    }
}

int huge_tlb_get_pagesize()
{
    HINSTANCE hDll;
    pfn_GetLargePageMinimum pGetLargePageMinimum;
    int pagesize = 0;

    // call succeeds only on Windows Server 2003 SP1 or later
    hDll = ::LoadLibrary(_T("kernel32.dll"));
    if (hDll == NULL) {
        huge_tlb_display_error(_T("LoadLibrary"), ::GetLastError());
        return -1;
    }

    pGetLargePageMinimum = (pfn_GetLargePageMinimum)::GetProcAddress(hDll, "GetLargePageMinimum");
    if (pGetLargePageMinimum == NULL) {
        huge_tlb_display_error(_T("GetProcAddress"), ::GetLastError());
        if (hDll)
            ::FreeLibrary(hDll);
        return -2;
    }

    pagesize = (*pGetLargePageMinimum)();

    if (hDll)
        ::FreeLibrary(hDll);
    return pagesize;
}

bool huge_tlb_adjust_token_privilege(TCHAR *pszPrivilege, BOOL bEnable)
{
    HANDLE              hToken = NULL;
    TOKEN_PRIVILEGES    tp = { 0 };
    BOOL                status;
    DWORD               error;
    bool                success = true;

    if (pszPrivilege == NULL) {
        huge_tlb_display_error(_T("pszPrivilege is NULL."), 0);
        return false;
    }

    // open process token
    if (!::OpenProcessToken(::GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, &hToken)) {
        huge_tlb_display_error(_T("OpenProcessToken"), ::GetLastError());
        success = false;
        goto HUGE_TLB_CLOSE_TOKEN;
    }

    if (hToken == NULL) {
        huge_tlb_display_error(_T("OpenProcessToken"), ::GetLastError());
        success = false;
        goto HUGE_TLB_CLOSE_TOKEN;
    }

    // get the luid
    if (!::LookupPrivilegeValue(NULL, pszPrivilege, &tp.Privileges[0].Luid)) {
        huge_tlb_display_error(_T("LookupPrivilegeValue"), ::GetLastError());
        success = false;
        goto HUGE_TLB_CLOSE_TOKEN;
    }

    tp.PrivilegeCount = 1;

    // enable or disable privilege
    if (bEnable)
        tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
    else
        tp.Privileges[0].Attributes = 0;

    // enable or disable privilege
    status = ::AdjustTokenPrivileges(hToken, FALSE, (PTOKEN_PRIVILEGES)&tp,
        sizeof(TOKEN_PRIVILEGES), (PTOKEN_PRIVILEGES)NULL, NULL);

    // It is possible for AdjustTokenPrivileges to return TRUE and still not succeed.
    // So always check for the last error value.
    error = ::GetLastError();
    if (!status || (error != ERROR_SUCCESS)) {
        huge_tlb_display_error(_T("AdjustTokenPrivileges"), ::GetLastError());
        success = false;
    }

HUGE_TLB_CLOSE_TOKEN:
    // close the handle
    if (hToken != NULL) {
        if (!CloseHandle(hToken)) {
            huge_tlb_display_error(_T("CloseHandle"), ::GetLastError());
            success = false;
        }
        else {
            hToken = NULL;
            success = true;
        }
    }
    return success;
}

/************************************************************************************
  Some remarks on VirtualAlloc and MEM_LARGE_PAGES

  From: http://blogs.msdn.com/b/oldnewthing/archive/2011/01/28/10121300.aspx

  Note that the MEM_LARGE_PAGES flag triggers an exception to the general principle
  that MEM_RESERVE only reserves address space, MEM_COMMIT makes the memory manager
  guarantee that physical pages will be there when you need them, and that
  the physical pages aren't actually allocated until you access the memory.
  Since very large pages have special physical memory requirements, the physical
  allocation is done up front so that the memory manager knows that when it comes
  time to produce the memory on demand, it can actually do so.
 ************************************************************************************/

void *huge_tlb_malloc(size_t size, size_t *real_alloc_size)
{
    return huge_tlb_malloc_ex(NULL, size, real_alloc_size);
}

#ifndef _MSC_VER

#define MEM_LARGE_PAGES  0x20000000
#define MEM_4MB_PAGES    0x80000000

#endif

void *huge_tlb_malloc_ex(void *base_address, size_t size, size_t *real_alloc_size)
{
    bool success;
    void *map_address = NULL;
    size_t alloc_size = 0;
    int large_pagesize;
    int huge_tlb_inited = s_huge_tlb_inited;
    if (huge_tlb_inited == 0)
        huge_tlb_inited = huge_tlb_init();
    if (huge_tlb_inited != 0) {
        large_pagesize = s_large_pagesize;
        if (size <= 0)
            alloc_size = large_pagesize;
        else
            alloc_size = ((size - 1) / large_pagesize + 1) * large_pagesize;
        success = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, TRUE);
        map_address = (void *)::VirtualAlloc(
            base_address,
            alloc_size,
            MEM_COMMIT | MEM_RESERVE | MEM_LARGE_PAGES,
            PAGE_READWRITE);
        if (map_address == NULL) {
            huge_tlb_display_error(_T("VirtualAlloc"), ::GetLastError());
            if (real_alloc_size != NULL)
                *real_alloc_size = 0;
            success = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, FALSE);
            return NULL;
        }
        else
            success = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, FALSE);
    }
    if (real_alloc_size != NULL)
        *real_alloc_size = alloc_size;
    return map_address;
}

bool huge_tlb_free(void *p)
{
    bool release_ok = false;
    if (p != NULL) {
        BOOL success;
        bool adjust_ok = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, TRUE);
        success = ::VirtualFree(p, 0, MEM_RELEASE);
        if (!success)
            huge_tlb_display_error(_T("VirtualFree"), ::GetLastError());
        else
            release_ok = true;
        adjust_ok = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, FALSE);
    }
    return release_ok;
}

bool huge_tlb_free_ex(void *p, size_t size)
{
    bool release_ok = false;
    if (p != NULL) {
        BOOL success;
        bool adjust_ok = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, TRUE);
        if (size == 0)
            success = ::VirtualFree(p, 0, MEM_RELEASE);
        else
            success = ::VirtualFree(p, size, /* MEM_LARGE_PAGES | */ MEM_DECOMMIT);
        if (!success) {
            huge_tlb_display_error(_T("VirtualFree"), ::GetLastError());
            printf("ptr = 0x%08X, size = %d\n\n", p, size);
            adjust_ok = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, FALSE);
        }
        else {
            adjust_ok = huge_tlb_adjust_token_privilege(SE_LOCK_MEMORY_NAME, FALSE);
            release_ok = true;
        }
    }
    return release_ok;
}
