#include <stdio.h>
#include <stdlib.h>
#include "tools.h"


/**
 * @brief Timing function that works for both Linux and MAC OS
 * 
 * @param ts Timespec variable.
 */
void my_gettime(struct timespec *ts) {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_MONOTONIC, ts);
#endif
}


/**
 * @brief Find elapsed time in seconds for given start timespec and
 * end timespec.
 * 
 * @param t_start Start timespec. 
 * @param t_end End timespec.
 * @return double 
 */
double elapsed_time(struct timespec *t_start, struct timespec *t_end)
{
    double time_secs = (t_end->tv_sec - t_start->tv_sec)
                     + (double) (t_end->tv_nsec - t_start->tv_nsec) * 1e-9;
    return time_secs;
}

