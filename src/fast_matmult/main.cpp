
// whether use crtdbg check?
#ifdef _DEBUG
#ifndef USE_CRTDBG_CHECK
#define USE_CRTDBG_CHECK    1
#endif
#endif

//
// From: http://msdn.microsoft.com/zh-cn/library/e5ewb1h3%28v=vs.90%29.aspx
// From: http://msdn.microsoft.com/en-us/library/x98tx3cf.aspx
//
#if USE_CRTDBG_CHECK
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#if _WIN32 || _WIN64
#include <tchar.h>
#include <objbase.h>
#include <mmsystem.h>
#endif
#include <malloc.h>

#if USE_CRTDBG_CHECK
#ifdef _DEBUG
// crtdbg.h must be behind the stdlib.h
#include <crtdbg.h>
#endif
#endif

#ifdef _DEBUG
#ifndef DEBUG_CLIENTBLOCK
#define DEBUG_CLIENTBLOCK   new(_CLIENT_BLOCK, __FILE__, __LINE__)
#endif
#define new                 DEBUG_CLIENTBLOCK
#else
#undef  DEBUG_CLIENTBLOCK
#define DEBUG_CLIENTBLOCK
#endif

#include "fast_matmult/fast_matmult.h"
#include "fast_matmult/huge_tlb.h"

using namespace annlab;

/**********************************************************
   Use Visual Leak Detector(vld) for Visual C++,
   Homepage: http://vld.codeplex.com/
 **********************************************************/
#ifdef _MSC_VER
#ifdef _DEBUG
// 如果你没有安装vld(Visual Leak Detector), 请注释掉这一句.
#include <vld.h>
#endif  /* _DEBUG */
#endif  /* _MSC_VER */

int main(int argc, char *argv[])
{
    unsigned int M, N, K;
    M = 256; N = 256; K = 256;
    //M = 512; N = 512; K = 512;
    M = 1024; N = 1024; K = 1024;
    //M = 2048; N = 2048; K = 2048;
    //M = 512; N = 1024; K = 256;

    printf("\n");
    HANDLE hCurrentProcess = GetCurrentProcess();
    DWORD dwProcessAffinity = 0, dwSystemAffinity = 0;
    DWORD dwAffinityMask = SetCPUAffinityMask4(0, 0, 0, 1);
    BOOL bAffResult;
    bAffResult = GetProcessAffinityMask(hCurrentProcess, &dwProcessAffinity, &dwSystemAffinity);
    if (bAffResult) {
        if (dwProcessAffinity != dwSystemAffinity)
            printf("This process can not utilize all processors.\n");

        while ((dwAffinityMask != 0) && (dwAffinityMask <= dwProcessAffinity)) {
            // Check to make sure we can utilize this processsor first.
            if ((dwAffinityMask & dwProcessAffinity) != 0) {
                bAffResult = SetProcessAffinityMask(hCurrentProcess, dwAffinityMask);
                if (bAffResult) {
                    Sleep(0);
                    //
                    printf("SetProcessAffinityMask(): dwAffinityMask = 0x%08X\n", dwAffinityMask);
                    DWORD dwProcessAffinityNew = 0;
                    bAffResult = GetProcessAffinityMask(hCurrentProcess, &dwProcessAffinityNew, &dwSystemAffinity);
                    if (dwProcessAffinityNew == dwAffinityMask) {
                        printf("SetProcessAffinityMask(): Success.\n");
                        bAffResult = SetThreadAffinityMask(GetCurrentThread(), dwAffinityMask);
                        break;
                    }
                }
            }
        }
    }
    printf("\n");

#if USE_CRTDBG_CHECK
    // 设置 CRT 报告模式
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);

    // 如果已经引用vld.h, 则不显示crtdbg的内存泄漏信息
#ifndef VLD_RPTHOOK_INSTALL
    // 进程退出时, 显示内存泄漏信息
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif  /* VLD_RPTHOOK_INSTALL */
#endif  /* USE_CRTDBG_CHECK */

    printf("Hello World!\n\n");

    int large_pagesize = huge_tlb_get_pagesize();
    int huge_tlb_inited = huge_tlb_init();
    if (huge_tlb_is_inited()) {
        printf("huge_tlb_init() is ok.\n\n");
    }

#if defined(USE_LARGE_PAGES) && (USE_LARGE_PAGES != 0)
#if 0
    TCHAR szName[] = _T("HUGE_TLB_LARGEPAGE");
    size_t alloc_size = 1024 * 1024 * sizeof(float_t) * 6;
    void *map_address = NULL;
  #if 0
    HANDLE hMapFile;
    hMapFile = ::CreateFileMapping(
        INVALID_HANDLE_VALUE,       // use paging file
        NULL,                       // default security
        PAGE_READWRITE | SEC_COMMIT | SEC_LARGE_PAGES,
        0,                          // max. object size
        alloc_size,                 // buffer size
        szName);                    // name of mapping object

    if (hMapFile == NULL)
        printf("CreateFileMapping(): error = %d\n\n", ::GetLastError());
    else
        printf("File mapping object successfulyl created.\n\n");

    map_address = (void *)::MapViewOfFile(
        hMapFile,                   // handle to map object
        FILE_MAP_ALL_ACCESS,        // read/write permission
        0,
        0,
        alloc_size);

    if (map_address) {
        printf("map_address = 0x%08X, alloc_size = 0x%08X (%d) byte(s).\n\n", map_address, alloc_size, alloc_size);
    }

    if (map_address == NULL)
        printf("MapViewOfFile(): error = %d\n\n", ::GetLastError());
    else
        printf("View of file successfully mapped.\n\n");

    system("pause");

    // do nothing, clean up an exit
    if (map_address)
        ::UnmapViewOfFile(map_address);
    if (hMapFile)
        ::CloseHandle(hMapFile);
  #else
    map_address = huge_tlb_malloc(alloc_size, &alloc_size);
    if (map_address) {
        printf("map_address = 0x%08X, alloc_size = 0x%08X (%d) byte(s).\n\n", map_address, alloc_size, alloc_size);
    }
    system("pause");

    if (map_address) {
        if (huge_tlb_free(map_address))
            printf("huge_tlb_free() done.\n\n");
        else
            printf("huge_tlb_free() fail.\n\n");
        map_address = NULL;
    }
  #endif
#endif
#endif  /* USE_LARGE_PAGES */

    matrix_matmult_test(M, K, N);

    huge_tlb_exit();

    printf("\n");
    system("pause");
    return 0;
}
