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

#define PT(l1, m1) ((m1) +((l1) *((l1) +1))/2)
#define YR(l2 , m2 ) (( m2 ) +( l2 ) +(( l2 ) *( l2) ) )
#define max(x,y) ((x) > (y) ? (x) : (y))

// TO PRINT IN COLOR
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

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


/**
 * @brief   Create a random matrix with each entry a random number within given range. 
 * 
 *          Note that each process within comm will have different random entries.
 *
 * @param Mat   Local part of the matrix.
 * @param m     Number of rows of the local copy of the matrix.
 * @param n     Number of columns of the local part of the matrix.
 */
void SetRandMat(double *Mat, int m, int n, double rand_min, double rand_max);


/**
 * @brief   Create a random matrix with each entry a random number within given range. 
 * 
 *          Note that each process within comm will have different random entries.
 *
 * @param Mat   Local part of the matrix.
 * @param m     Number of rows of the local copy of the matrix.
 * @param n     Number of columns of the local part of the matrix.
 */
void SetRandMat_complex(double complex *Mat, int m, int n, double rand_min, double rand_max);


/**
 * @brief Compare two double arrays and return the maximum difference.
 * 
 * @param a First array.
 * @param b Second array.
 * @param n Length of the array a and array b.
 * @param tol Tolerance for maximum difference.
 * @return double 
 */
double check_double_arrays(
    const double *a, const int stride_a, const double *b, const int stride_b,
    const int n, const double tol);


/**
 * @brief Convert Cartesian coordiantes to polar coordinates.
 * 
 * @param n Number of coordinates to convert.
 * @param x X coordinates.
 * @param y Y coordinates.
 * @param z Z coordinates.
 * @param r Radius.
 * @param theta Angle between z axis and the position vector.
 * @param phi Angle between x axis and the projection of position vector on x-y plane.
 */
void cart2polar(const int n, const double *x, const double *y, const double *z,
	double *r, double *theta, double *phi);


void show_results(const int n, const double *x, const double *y, const double *z,
	const int Lmax, const double *Ylm, const double *Y);


#endif // _TOOLS_H
