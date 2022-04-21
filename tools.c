#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <assert.h>
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


/**
 * @brief   Create a random matrix with each entry a random number within given range. 
 * 
 *          Note that each process within comm will have different random entries.
 *
 * @param Mat   Local part of the matrix.
 * @param m     Number of rows of the local copy of the matrix.
 * @param n     Number of columns of the local part of the matrix.
 */
void SetRandMat(double *Mat, int m, int n, double rand_min, double rand_max)
{
	int rank = 0, i, len_tot;
	int seed_shift = 1;
	int seed_temp = rank * 100 + seed_shift;
	srand(seed_temp);

	len_tot = m * n;
	for (i = 0; i < len_tot; i++) {
		Mat[i] = rand_min + (rand_max - rand_min) * ((double) rand() / RAND_MAX);
	}
}


/**
 * @brief   Create a random matrix with each entry a random number within given range. 
 * 
 *          Note that each process within comm will have different random entries.
 *
 * @param Mat   Local part of the matrix.
 * @param m     Number of rows of the local copy of the matrix.
 * @param n     Number of columns of the local part of the matrix.
 */
void SetRandMat_complex(double complex *Mat, int m, int n, double rand_min, double rand_max)
{
	int rank = 0, i, len_tot;
	int seed_shift = 1;
	int seed_temp = rank * 100 + seed_shift;
	srand(seed_temp);

	len_tot = m * n;
	for (i = 0; i < len_tot; i++) {
		Mat[i] = rand_min + (rand_max - rand_min) * ((double) rand() / RAND_MAX);
	}
}


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
	const int n, const double tol)
{
	// double tol = 1e-8;
	double err = 0.0;
	for (int i = 0; i < n; i++) {
		err = max(err, fabs(a[i*stride_a] - b[i*stride_b]));
	}
	// if (err >= tol)
	//     printf("In function check_double_arrays: err = %.3e\n",err);
	// return (int) (err >= tol); // 1 - error, 0 - success
	return err;
}


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
	double *r, double *theta, double *phi)
{
	for (int i = 0; i < n; i++) {
		r[i] = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
		theta[i] = atan(sqrt(x[i]*x[i] + y[i]*y[i])/z[i]);
		if (theta[i]<0.0) theta[i] += M_PI;
		if (z[i]==0) theta[i] = M_PI/2.0;
		phi[i] = atan(y[i]/x[i]);
		if (x[i] < 0.0) phi[i] += M_PI;  // 2nd Quandrant
		if (x[i] > 0.0 && y[i] < 0.0) phi[i] += 2*M_PI;  // 4th Quadrant
		if (x[i] == 0.0 && y[i] < 0.0) phi[i] = 3*M_PI/2;
		if (x[i]==0.0 && y[i] >0.0) phi[i] = M_PI/2.0;
		// if (x[i]>0.0 && y[i]>0.0) phi[i] = phi[i];  // 1st Quandrant
		// if (x[i]<0.0 && y[i]>0.0) phi[i] = M_PI+phi[i];  // 2nd Quandrant
		// if (x[i]<0.0 && y[i]<0.0) phi[i] = M_PI+phi[i];  // 3rd Quadrant
		// if (x[i]>0.0 && y[i]<0.0) phi[i] = 2.0*M_PI+phi[i];  // 4th Quadrant
		// if (x[i]==0.0 && y[i] >0.0) phi[i] = M_PI/2.0;
		// if (x[i]==0.0 && y[i] <0.0) phi[i] = 3.0*M_PI/2.0;
	}
}


/**
 * @brief Show the results that are not correct (error > tol).
 * 
 * @param n Number of coordinates.
 * @param x X coords.
 * @param y Y coords.
 * @param z Z coords.
 * @param Lmax 
 * @param Ylm 
 * @param Y 
 */
void show_results(const int n, const double *x, const double *y, const double *z,
	const int Lmax, const double *Ylm, const double *Y)
{
	int memsize_Y = (Lmax+1)*(Lmax+1);
	// double *r = (double *) malloc(n * sizeof(*r));
	// double *theta = (double *) malloc (n * sizeof(*theta));
	// double *phi = (double *) malloc (n * sizeof(*phi));
	// assert(r != NULL && theta != NULL && phi != NULL);
	// cart2polar(n, x, y, z, r, theta, phi);
	// free(r); free(theta); free(phi);

	double error = 0.0, tol = 1e-8;
	for (int l = 0; l < Lmax; l++) {
		for (int m = -l; m <= l; m++) {
			// double error_lm = check_double_arrays(Ylm, 1, &Y[PT(l, m)], memsize_Y, n, tol);
			// double error_lm = check_double_arrays(Ylm, 1, &Y[YR(l, m)], memsize_Y, n, tol);
			double error_lm = check_double_arrays(Ylm, 1, &Y[YR(l, m)], memsize_Y, n, tol);
			error = max(error, error_lm);
			if (error_lm > tol) {
				for (int i = 0; i < n; i++) {
					printf("(x,y,z) = (%.3f,%.3f,%.3f), l = %2d, m = %2d, Ylm = %.3f, Y[l,m] = %.3f, error_lm = %.3e\n",
						x[i],y[i],z[i],l,m,Ylm[i],Y[YR(l, m)+i*memsize_Y],error_lm);
				}
			}
		}
	}
}

