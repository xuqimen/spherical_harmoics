/**
 * @file test.c
 * @author Qimen Xu (qimenxu@foxmail.com)
 *         Shashikant Kumar (shashikant@gatech.edu)
 * @brief Tests for spherical harmonics routines.
 * @version 0.1
 * @date 2022-04-19
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "tools.h"
#include "spherical_harmonics.h"

int main() {
	double theta, phi;
	double complex *Y, *dY_th, *dY_phi;
	int Lmax, memsize_Y;
    struct timespec t_start, t_end, t1, t2;
    double time_secs;

    my_gettime(&t1);
 
	Lmax=2;
	theta=1.2;
	phi = 2.1;
	memsize_Y = (Lmax+1)*(Lmax+1);
	
    my_gettime(&t_start);
    Y = (double complex *) malloc(sizeof(double complex)*memsize_Y);
	dY_phi = (double complex *) malloc(sizeof(double complex)*memsize_Y);
	dY_th = (double complex *) malloc(sizeof(double complex)*memsize_Y);
    my_gettime(&t_end);
    time_secs = elapsed_time(&t_start, &t_end);
    printf("Run-time of the memory allocation: %.3lf ms\n", time_secs*1000.0);
 
    my_gettime(&t_start);
    sph_harmonics(theta, phi, Lmax, Y, dY_th, dY_phi);
    my_gettime(&t_end);
    // time in seconds
    time_secs = elapsed_time(&t_start, &t_end);
    // double time_secs = (t_end.tv_sec - t_start.tv_sec)
    //                 + (double) (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    printf("Run-time of the sph_harmonics: %.3lf ms\n", time_secs*1000.0);
    // fprintf(stderr, "Run-time of the sph_harmonics: %8.0lf ms\n", time_secs*1000.0);
 
	// for (int l=0; l <= Lmax; l++){
	// 	for (int m=-l; m<=l; m++){
	// 		printf("Y(%d,%d,%f,%f): %f + %f i,\n dY_theta(%d,%d,%f,%f): %f + %f i,\n dY_phi(%d,%d,%f,%f): %f + %f i\n",
	// 			  l,m,theta,phi, creal(Y[YR(l,m)]), cimag(Y[YR(l,m)]),  l,m,theta,phi, creal(dY_th[YR(l,m)]), cimag(dY_th[YR(l,m)]),
	// 			  l,m,theta,phi, creal(dY_phi[YR(l,m)]), cimag(dY_phi[YR(l,m)]));
	// 			  printf("\n");
	// 	}
	// }
 
	free(Y);
	free(dY_phi);
	free(dY_th);

    my_gettime(&t2);
    time_secs = elapsed_time(&t1, &t2);
    printf("Run-time of total program: %.3lf ms\n", time_secs*1000.0);

	return 0;
}
