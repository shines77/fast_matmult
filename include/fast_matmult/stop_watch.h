/*
    Copyright 2005-2011 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#ifndef _ANNLAB_STOP_WATCH_H_
#define _ANNLAB_STOP_WATCH_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#ifdef _MSC_VER
#include <fast_matmult/vs/stdint.h>
#else
#include <stdint.h>
#endif

#if _WIN32 || _WIN64
#include <windows.h>
#include <time.h>
#elif __linux__
#include <ctime>
#else /* generic Unix */
#include <sys/time.h>
#endif /* (choice of OS) */

namespace annlab {

//! Absolute timestamp
/** @ingroup timing */
class stop_watch
{
public:
    typedef int64_t timestamp_t;

public:
    //! Construct an absolute timestamp initialized to zero.
    stop_watch() : bIsRunning(false), startTime(0), stopTime(0), elapsedTime(0), elapsedTimeTotal(0) {};
    stop_watch(const stop_watch &src);

    stop_watch &operator =(const stop_watch &t);

    bool    isRunning(void);

    void    clear(void);
    void    reset(void);

    //! restart() is equivalent to reset() and begin()
    void    restart(void);
    void    reset_and_begin(void);

    void    start(void);
    void    stop(void);

    void    begin(void);
    void    end(void);

    //! Return current time.
    static timestamp_t  now(void);

    //! Return current time(double).
    static double       nowf(void);

    static double       intervalSeconds(int64_t t1, int64_t t2);
    static double       intervalSeconds(double t1, double t2);

    //! Return current time(millisecs).
    static timestamp_t  currentTimeMillis(void);
    static double       currentTimeMillisf(void);

    double  getSeconds(void);
    double  getMillisec(void);

    double  getTotalSeconds(void);
    double  getTotalMillisec(void);

    double  getUsedTime(void);
    double  getUsedTimeTotal(void);

protected:
    void    native_start();
    void    native_stop();

    static timestamp_t  native_now();
    static double       native_nowf();

private:
    bool        bIsRunning;
    timestamp_t startTime, stopTime;
    timestamp_t elapsedTime;
    timestamp_t elapsedTimeTotal;
};

inline stop_watch::stop_watch(const stop_watch &src)
{
    bIsRunning          = src.bIsRunning;
    startTime           = src.startTime;
    stopTime            = src.stopTime;
    elapsedTime         = src.elapsedTime;
    elapsedTimeTotal    = src.elapsedTimeTotal;
}

///*
inline stop_watch &stop_watch::operator =(const stop_watch &t)
{
    bIsRunning          = t.bIsRunning;
    startTime           = t.startTime;
    stopTime            = t.stopTime;
    elapsedTime         = t.elapsedTime;
    elapsedTimeTotal    = t.elapsedTimeTotal;
    return *this;
}
//*/

inline void stop_watch::native_start()
{
#if _WIN32 || _WIN64
    LARGE_INTEGER qp_cnt;
    QueryPerformanceCounter(&qp_cnt);
    startTime = static_cast<timestamp_t>(qp_cnt.QuadPart);
    stopTime  = startTime;
#elif __linux__
    struct timespec ts;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        clock_gettime(CLOCK_REALTIME, &ts);
    _GMTL_ASSERT(status == 0, "CLOCK_REALTIME not supported");
    startTime = static_cast<timestamp_t>(static_cast<int64_t>(1000000000UL) * static_cast<int64_t>(ts.tv_sec) + static_cast<int64_t>(ts.tv_nsec));
    stopTime  = startTime;
#else /* generic Unix */
    struct timeval tv;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        gettimeofday(&tv, NULL);
    _GMTL_ASSERT(status == 0, "gettimeofday failed");
    startTime = static_cast<timestamp_t>(static_cast<int64_t>(1000000) * static_cast<int64_t>(tv.tv_sec) + static_cast<int64_t>(tv.tv_usec));
    stopTime  = startTime;
#endif /*(choice of OS) */
}

inline void stop_watch::native_stop()
{
#if _WIN32 || _WIN64
    LARGE_INTEGER qp_cnt;
    QueryPerformanceCounter(&qp_cnt);
    stopTime = static_cast<timestamp_t>(qp_cnt.QuadPart);
#elif __linux__
    struct timespec ts;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        clock_gettime(CLOCK_REALTIME, &ts);
    _GMTL_ASSERT(status == 0, "CLOCK_REALTIME not supported");
    stopTime = static_cast<timestamp_t>(static_cast<int64_t>(1000000000UL) * static_cast<int64_t>(ts.tv_sec) + static_cast<int64_t>(ts.tv_nsec));
#else /* generic Unix */
    struct timeval tv;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        gettimeofday(&tv, NULL);
    _GMTL_ASSERT(status == 0, "gettimeofday failed");
    stopTime = static_cast<timestamp_t>(static_cast<int64_t>(1000000) * static_cast<int64_t>(tv.tv_sec) + static_cast<int64_t>(tv.tv_usec));
#endif /*(choice of OS) */
}

inline stop_watch::timestamp_t stop_watch::native_now(void)
{
    timestamp_t result;

#if _WIN32 || _WIN64
    LARGE_INTEGER qp_cnt, qp_freq;
    QueryPerformanceCounter(&qp_cnt);
    QueryPerformanceFrequency(&qp_freq);
    result = static_cast<timestamp_t>(((double)qp_cnt.QuadPart / (double)qp_freq.QuadPart) * 1000000000.0);
#elif __linux__
    struct timespec ts;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        clock_gettime(CLOCK_REALTIME, &ts);
    _GMTL_ASSERT(status == 0, "CLOCK_REALTIME not supported");
    result = static_cast<timestamp_t>(static_cast<int64_t>(1000000000UL) * static_cast<int64_t>(ts.tv_sec) + static_cast<int64_t>(ts.tv_nsec));
#else /* generic Unix */
    struct timeval tv;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        gettimeofday(&tv, NULL);
    _GMTL_ASSERT(status == 0, "gettimeofday failed");
    result = static_cast<timestamp_t>(static_cast<int64_t>(1000000000UL) * static_cast<int64_t>(tv.tv_sec) + static_cast<int64_t>(1000UL) * static_cast<int64_t>(tv.tv_usec));
#endif /*(choice of OS) */

    return result;
}

inline double stop_watch::native_nowf(void)
{
    double result;

#if _WIN32 || _WIN64
    LARGE_INTEGER qp_cnt, qp_freq;
    QueryPerformanceCounter(&qp_cnt);
    QueryPerformanceFrequency(&qp_freq);
    result = (double)qp_cnt.QuadPart / (double)qp_freq.QuadPart;
#elif __linux__
    int64_t time_usecs;
    struct timespec ts;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        clock_gettime(CLOCK_REALTIME, &ts);
    _GMTL_ASSERT(status == 0, "CLOCK_REALTIME not supported");
    time_usecs = static_cast<int64_t>(1000000000UL) * static_cast<int64_t>(ts.tv_sec) + static_cast<int64_t>(ts.tv_nsec);
    result = (double)time_usecs * 1E-9
#else /* generic Unix */
    int64_t time_usecs;
    struct timeval tv;
#if GMTL_USE_ASSERT
    int status =
#endif /* GMTL_USE_ASSERT */
        gettimeofday(&tv, NULL);
    _GMTL_ASSERT(status == 0, "gettimeofday failed");
    time_usecs = static_cast<int64_t>(1000000UL) * static_cast<int64_t>(tv.tv_sec) + static_cast<int64_t>(tv.tv_usec);
    result = (double)time_usecs * 1E-6;
#endif /*(choice of OS) */

    return result;
}

inline bool stop_watch::isRunning(void)
{
    return bIsRunning;
}

inline void stop_watch::clear(void)
{
    elapsedTime = 0;
    elapsedTimeTotal = 0;

    native_start();

    bIsRunning = false;
}

inline void stop_watch::reset(void)
{
    elapsedTime = 0;

    native_start();

    bIsRunning = false;
}

//! restart() is equivalent to reset() and begin()
inline void stop_watch::restart(void)
{
    elapsedTime = 0;

    native_start();

    bIsRunning = true;
}

inline void stop_watch::reset_and_begin(void)
{
    restart();
}

inline void stop_watch::start(void)
{
    native_start();

    bIsRunning = true;
}

inline void stop_watch::begin(void)
{
    start();
}

inline void stop_watch::stop(void)
{
    native_stop();

    if (bIsRunning) {
        elapsedTime = stopTime - startTime;
        if (elapsedTime >= 0)
            elapsedTimeTotal += elapsedTime;
        else
            elapsedTimeTotal -= elapsedTime;
    }

    bIsRunning = false;
}

inline void stop_watch::end(void)
{
    stop();
}

inline stop_watch::timestamp_t stop_watch::now(void)
{
    return native_now();
}

inline double stop_watch::nowf(void)
{
    return native_nowf();
}

inline double stop_watch::intervalSeconds(int64_t t1, int64_t t2)
{
    double seconds = (double)(t2 - t1) * 1E-9;
    return seconds;
}

inline double stop_watch::intervalSeconds(double t1, double t2)
{
    double seconds = (double)(t2 - t1);
    return seconds;
}

inline stop_watch::timestamp_t stop_watch::currentTimeMillis(void)
{
    timestamp_t now_usecs = stop_watch::native_now();
    return now_usecs / static_cast<timestamp_t>(1000UL);
}

inline double stop_watch::currentTimeMillisf(void)
{
    double now_usecs = stop_watch::native_nowf();
    return now_usecs * 1E-3;
}

inline double stop_watch::getSeconds(void)
{
    if (bIsRunning)
        stop();

#if _WIN32 || _WIN64
    LARGE_INTEGER qp_freq;
    QueryPerformanceFrequency(&qp_freq);
    return (double)elapsedTime / (double)qp_freq.QuadPart;
#elif __linux__
    return (double)elapsedTime * 1E-9;
#else /* generic Unix */
    return (double)elapsedTime * 1E-6;
#endif /* (choice of OS) */
}

inline double stop_watch::getMillisec(void)
{
    if (bIsRunning)
        stop();

#if _WIN32 || _WIN64
    LARGE_INTEGER qp_freq;
    QueryPerformanceFrequency(&qp_freq);
    return ((double)elapsedTime / (double)qp_freq.QuadPart) * 1E3;
#elif __linux__
    return (double)elapsedTime * 1E-6;
#else /* generic Unix */
    return (double)elapsedTime * 1E-3;
#endif /* (choice of OS) */
}

inline double stop_watch::getTotalSeconds(void)
{
    if (bIsRunning)
        stop();

#if _WIN32 || _WIN64
    LARGE_INTEGER qp_freq;
    QueryPerformanceFrequency(&qp_freq);
    return (double)elapsedTimeTotal / (double)qp_freq.QuadPart;
#elif __linux__
    return (double)elapsedTimeTotal * 1E-9;
#else /* generic Unix */
    return (double)elapsedTimeTotal * 1E-6;
#endif /* (choice of OS) */
}

inline double stop_watch::getTotalMillisec(void)
{
    if (bIsRunning)
        stop();

#if _WIN32 || _WIN64
    LARGE_INTEGER qp_freq;
    QueryPerformanceFrequency(&qp_freq);
    return ((double)elapsedTimeTotal / (double)qp_freq.QuadPart) * 1E3;
#elif __linux__
    return (double)elapsedTimeTotal * 1E-6;
#else /* generic Unix */
    return (double)elapsedTimeTotal * 1E-3;
#endif /* (choice of OS) */
}

inline double stop_watch::getUsedTime(void)
{
    return getSeconds();
}

inline double stop_watch::getUsedTimeTotal(void)
{
    return getTotalSeconds();
}

}  /* namespace annlab */

#endif  /* _ANNLAB_STOP_WATCH_H_ */
