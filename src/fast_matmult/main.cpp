
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define _N   2048
double A[_N][_N];
double B[_N][_N];
double C[_N][_N];

void simple_martrix_test()
{
    // populate the matrices with random values between 0.0 and 1.0
    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            A[i][j] = (double)rand() / (double)RAND_MAX;
            B[i][j] = (double)rand() / (double)RAND_MAX;
            C[i][j] = 0.0;
        }
    }

    stop_watch sw;
    double time_spent;

    // matrix multiplication
    sw.start();
    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            for (int k = 0; k < _N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    sw.stop();

    time_spent = sw.getMillisec();
    printf("C[1024][1024] = %0.5f\n\n"
           "Elapsed time: %6.2f ms\n\n",
           C[1024][1024], time_spent);
}

/**********************************************************
 *
 *  Use Visual Leak Detector(vld) for Visual C++,
 *  Homepage: http://vld.codeplex.com/
 *
 **********************************************************/
#ifdef _MSC_VER
#ifdef _DEBUG
// �����û�а�װvld(Visual Leak Detector), ��ע�͵���һ��.
//#include <vld.h>
#endif  /* _DEBUG */
#endif  /* _MSC_VER */

/* ����CRTDBG�Ļ���(Debugģʽ��, ����ڴ�Խ����ڴ�й©����) */

void set_crtdbg_env()
{
/* ʹ��CRTDBG�������ڴ�Խ������, �����ʹ����vld, �ڴ�й©��Ϣ�ɹر� */
#if USE_CRTDBG_CHECK
    // ���� CRT ����ģʽ
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);

    // ����Ѿ�����vld.h, ����ʾcrtdbg���ڴ�й©��Ϣ
#ifndef VLD_RPTHOOK_INSTALL
    // �����˳�ʱ, ��ʾ�ڴ�й©��Ϣ
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
    tips_format_text[1] = "����ѡ����: [�˳� = %d]: ? ";

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
                // ���������ַ�����ָ���ķ�Χ, ��ʾһ��(500����)���ָ̻����һ����ȷ��ѡ��ֵ
                printf("\x08%d", tmp_value);
                fflush(stdout);

                // ����500����
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
        "��ѡ����Ҫ���е�ģʽ:\n\n"
        ""
        "[1] = ��ʹ�� tiling �ֿ鼼���Ĵ� C/C++ ����.\n"
        "[2] =   ʹ�� tiling �ֿ鼼���Ĵ� C/C++ ����.\n"
        "[3] =   ʹ�� tiling �ֿ鼼����ʹ�� SSEx ָ���Ż�.\n"
        "\n"
        "[4] = ���� �� C/C++ ����� ����: ���� [1] [2].\n"
        "[5] = ���� ʹ�÷ֿ鼼����  ����: ���� [2] [3].\n"
        "[6] = ���� �����оٵ�      ����: ���� [1] [2] [3]. (Ĭ��)\n"
        "\n"
        "[0] = �˳�����.\n\n"
        ""
        "����������ѡ���Իس�������...\n\n";

    tips_format_text[1] = "����ѡ����: [%d = �˳�]: ? ";

    return get_user_choice(lang_id, display_text[lang_id], tips_format_text[lang_id], 0, 6, default_value);
}

/* Ԥ��ʱ������Ҫ����500����, ���������, ��������������СԤ��ʱ�� */
void iso_cpu_warm_up()
{
#ifndef _DEBUG
    stop_watch sw;
    volatile int sum = 0;
    double elapsedTime = 0.0;
    printf("CPU warm up start ...\n");
    do {
        sw.restart();
        // ����д����ı������ܷ�������һ���̶�ֵ��ţ����, Ӧ��û��
        for (int i = 0; i < 10000; ++i) {
            sum += i;
            // ѭ��˳�����ߵ�������
            for (int j = 5000; j >= 0; --j) {
                sum -= j;
            }
        }
        sw.stop();
        elapsedTime += sw.getMillisec();
    } while (elapsedTime < 500.0);
    // ���sum��ֵֻ��Ϊ�˷�ֹ��������ѭ���Ż���
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

    // ����CRTDBG�Ļ���(Debugģʽ��, ����ڴ�Խ����ڴ�й©����)
    set_crtdbg_env();

    // ���ý��̺��̵߳���Ե��
    set_thread_affinity();

    // ��ȡ�û�����ĳ�������ģʽroutine_mode
    //routine_mode = 3;
    routine_mode = get_routine_mode(3);
    if (routine_mode == GETCH_EXIT_PROGRAM) {
        system("pause");
        return 0;
        //goto L_EXIT_MAIN;
    }

    printf("\n");

    //simple_martrix_test();

#if 1

    if (lang_id == LANG_ZH_CN) {
        printf("����������ά��: ");
        printf("[0 = exit]\n");
        printf("(������2���ݴη�, ����: 1024, ���ֵΪ8192)\n\n");
    }
    else {
	    printf("Please enter the dimension of the matrix: ");
        printf("[0 = exit]\n");
        printf("(must be power of 2, for example: 1024, maximum value is 8192)\n\n");
    }

	printf("Dim = ? ");
#ifdef _MSC_VER
	scanf_s("%u", &n);
#else
    scanf("%u", &n);
#endif
    printf("\n");

	if (n == 0) {
        system("pause");
	    return 0;
		//goto L_EXIT_MAIN;
    }

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

//L_EXIT_MAIN:
    system("pause");
    return 0;
}
