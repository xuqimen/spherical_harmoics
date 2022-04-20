/**
 * @file tools.h
 * @author Qimen Xu (qimenxu@foxmail.com)
 * @brief Tool functions.
 * @version 0.1
 * @date 2022-04-19
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef _TOOLS_H
#define _TOOLS_H


#include <time.h>
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

/**
 * @brief Timing function that works for both Linux and MAC OS
 * 
 * @param ts Timespec variable.
 */
void my_gettime(struct timespec *ts);


/**
 * @brief Find elapsed time in seconds for given start timespec and
 * end timespec.
 * 
 * @param t_start Start timespec. 
 * @param t_end End timespec.
 * @return double 
 */
double elapsed_time(struct timespec *t_start, struct timespec *t_end);


#endif
