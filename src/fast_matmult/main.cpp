
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
//#include <tchar.h>
#include <objbase.h>
#include <mmsystem.h>
#endif
#include <conio.h>
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
#include "fast_matmult/get_char.h"
#include "fast_matmult/huge_tlb.h"
#include "fast_matmult/stop_watch.h"

using namespace annlab;

/**********************************************************
   Use Visual Leak Detector(vld) for Visual C++,
   Homepage: http://vld.codeplex.com/
 **********************************************************/
#ifdef _MSC_VER
#ifdef _DEBUG
// 如果你没有安装vld(Visual Leak Detector), 请注释掉这一句.
//#include <vld.h>
#endif  /* _DEBUG */
#endif  /* _MSC_VER */

void set_crtdbg_env()
{
/* 使用CRTDBG将会检查内存越界问题, 如果你使用了vld, 内存泄漏信息可关闭 */
#if USE_CRTDBG_CHECK
    // 设置 CRT 报告模式
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);

    // 如果已经引用vld.h, 则不显示crtdbg的内存泄漏信息
#ifndef VLD_RPTHOOK_INSTALL
    // 进程退出时, 显示内存泄漏信息
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif  /* VLD_RPTHOOK_INSTALL */
#endif  /* USE_CRTDBG_CHECK */
}

static size_t __cdecl _next_power_of_2(size_t x)
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

int get_user_choice(int lang_id, char *display_text, char *tips_format_text_,
                    int min_value, int max_value, int default_value)
{
    const int exit_value    = GETCH_EXIT_PROGRAM;
    int input_value         = GETCH_DEFUALT_VALUE;
    int tmp_value;

    char *tips_format_text[MAX_LANG_ID];

    tips_format_text[0] = "Your choice is: [exit = %d]: ? ";
    tips_format_text[1] = "您的选择是: [退出 = %d]: ? ";

    printf("%s", display_text);
    if (tips_format_text_ == NULL)
        printf(tips_format_text[lang_id], exit_value);
    else
        printf(tips_format_text_, exit_value);
    printf("%d", default_value);

    int nchar;
    do {
        nchar = iso_getch();
        if (nchar == MM_VT_KEY_RETURN) {
            if (input_value == GETCH_DEFUALT_VALUE)
                input_value = default_value;
            break;
        }
        else if (nchar >= '0' && nchar <= '9') {
            tmp_value = nchar - '0';
            if (tmp_value >= min_value && tmp_value <= max_value) {
                input_value = tmp_value;
                printf("\x08%d", input_value);
                fflush(stdout);
            }
            else {
                // 如果输入的字符不在指定的范围, 显示一下(500毫秒)立刻恢复最后一次正确的选项值
                printf("\x08%d", tmp_value);
                fflush(stdout);

                // 休眠500毫秒
                iso_sleep(500);

                if (input_value != GETCH_DEFUALT_VALUE)
                    printf("\x08%d", input_value);
                else
                    printf("\x08%d", default_value);
                fflush(stdout);
            }
        }
    } while (1);

    printf("\n\n");
    return input_value;
}

int get_routine_mode(int default_value)
{
    int lang_id = get_current_langid();

    char *display_text[MAX_LANG_ID];
    char *tips_format_text[MAX_LANG_ID];

    display_text[0] =
        "Please choice program's routine mode:\n\n"
        ""
        "[1] = Pure C/C++ code   that unused tiling.\n"
        "[2] = Pure C/C++ code   that    use tiling.\n"
        "[3] = Use SSEx instructions and use tiling.\n"
        "\n"
        "[4] = All of pure C/C++ codes:          [1] and [2].\n"
        "[5] = All of the tests that use tiling: [2] [3].\n"
        "[6] = All of the above tests:           [1] [2] [3]. (default)\n"
        "\n"
        "[0] = Exit program.\n\n"
        ""
        "Please input your choice and press enter key to continue...\n\n";

    tips_format_text[0] = "Your choice is: [%d = exit]: ? ";

    display_text[1] =
        "请选择你要运行的模式:\n\n"
        ""
        "[1] = 不使用 tiling 分块技术的纯 C/C++ 代码.\n"
        "[2] =   使用 tiling 分块技术的纯 C/C++ 代码.\n"
        "[3] =   使用 tiling 分块技术和使用 SSEx 指令优化.\n"
        "\n"
        "[4] = 所有 纯 C/C++ 代码的 测试: 包含 [1] [2].\n"
        "[5] = 所有 使用分块技术的  测试: 包含 [2] [3].\n"
        "[6] = 所有 以上列举的      测试: 包含 [1] [2] [3]. (默认)\n"
        "\n"
        "[0] = 退出程序.\n\n"
        ""
        "请输入您的选择并以回车键结束...\n\n";

    tips_format_text[1] = "您的选择是: [%d = 退出]: ? ";

    return get_user_choice(lang_id, display_text[lang_id], tips_format_text[lang_id], 0, 6, default_value);
}

/* 预热时间至少要大于500毫秒, 如果还不够, 可以自行增加最小预热时间 */
void iso_cpu_warm_up()
{
#ifndef _DEBUG
    stop_watch sw;
    volatile int sum = 0;
    double elapsedTime = 0.0;
    printf("CPU warm up start ...\n");
    do {
        sw.restart();
        // 如果有聪明的编译器能发现这是一个固定值就牛比了, 应该没有
        for (int i = 0; i < 10000; ++i) {
            sum += i;
            // 循环顺序故意颠倒过来的
            for (int j = 5000; j >= 0; --j) {
                sum -= j;
            }
        }
        sw.stop();
        elapsedTime += sw.getMillisec();
    } while (elapsedTime < 500.0);
    // 输出sum的值只是为了防止编译器把循环优化掉
    printf("sum = %u, time: %0.3f ms\n", sum, elapsedTime);
    printf("CPU warm up done ... \n\n");
#endif
}

void set_thread_affinity()
{
    bool echo = false;
#if _WIN32 || _WIN64
    if (echo)
        printf("\n");
    HANDLE hCurrentProcess = GetCurrentProcess();
    DWORD dwProcessAffinity = 0, dwSystemAffinity = 0;
    DWORD dwAffinityMask = SetCPUAffinityMask4(1, 0, 0, 0);
    BOOL bAffResult;
    bAffResult = GetProcessAffinityMask(hCurrentProcess, &dwProcessAffinity, &dwSystemAffinity);
    if (bAffResult) {
        if (dwProcessAffinity != dwSystemAffinity) {
            if (echo)
                printf("This process can not utilize all processors.\n");
        }

        while ((dwAffinityMask != 0) && (dwAffinityMask <= dwProcessAffinity)) {
            // Check to make sure we can utilize this processsor first.
            if ((dwAffinityMask & dwProcessAffinity) != 0) {
                bAffResult = SetProcessAffinityMask(hCurrentProcess, dwAffinityMask);
                if (bAffResult) {
                    // Wait for the process affinity effected
                    Sleep(0);
                    //
                    if (echo)
                        printf("SetProcessAffinityMask(): dwAffinityMask = 0x%08X\n", dwAffinityMask);
                    DWORD dwProcessAffinityNew = 0;
                    bAffResult = GetProcessAffinityMask(hCurrentProcess, &dwProcessAffinityNew, &dwSystemAffinity);
                    if (dwProcessAffinityNew == dwAffinityMask) {
                        if (echo)
                            printf("SetProcessAffinityMask(): Success.\n");
                        bAffResult = SetThreadAffinityMask(GetCurrentThread(), dwAffinityMask);
                        Sleep(0);
                        break;
                    }
                }
            }
        }
    }
    if (echo)
        printf("\n");
#endif
}

int main(int argc, char *argv[])
{
    unsigned int M, N, K, n = 1024;
    bool echo = false;
    int routine_mode = 0;

    int lang_id;
    lang_id = set_current_langid(LANG_USER);

    // 设置CRTDBG的环境(Debug模式下, 检查内存越界和内存泄漏问题)
    set_crtdbg_env();

    // 设置进程和线程的亲缘性
    set_thread_affinity();

    // 获取用户输入的程序运行模式routine_mode
    //routine_mode = 3;
    routine_mode = get_routine_mode(3);
    if (routine_mode == GETCH_EXIT_PROGRAM)
        goto _EXIT_MAIN;

#if 1

    if (lang_id == LANG_ZH_CN) {
        printf("请输入矩阵的维数: ");
        printf("[0 = exit]\n");
        printf("(必须是2的幂次方, 例如: 1024, 最大值为8192)\n\n");
    }
    else {
	    printf("Please enter the dimension of the matrix: ");
        printf("[0 = exit]\n");
        printf("(must be power of 2, for example: 1024, maximum value is 8192)\n\n");
    }

	printf("Dim = ? ");
	scanf_s("%u", &n);
    printf("\n");

	if (n == 0)
		goto _EXIT_MAIN;

    // n round to power of 2
    n = _next_power_of_2(n);

    if (n < 64)
        n = 64;

    if (n > 8192)
        n = 8192;

    M = n; N = n; K = n;
#else
    //M = 256; N = 256; K = 256;
    //M = 512; N = 512; K = 512;
    M = 1024; N = 1024; K = 1024;
    //M = 2048; N = 2048; K = 2048;
#endif

    printf("M = %d, N = %d, K = %d\n\n", M, N, K);

    int large_pagesize = huge_tlb_get_pagesize();
    int huge_tlb_inited = huge_tlb_init();
    if (huge_tlb_is_inited()) {
        if (echo)
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

    iso_cpu_warm_up();

    matrix_matmult_test(routine_mode, M, K, N);

    huge_tlb_exit(echo);
    printf("\n");

_EXIT_MAIN:
    system("pause");
    return 0;
}
