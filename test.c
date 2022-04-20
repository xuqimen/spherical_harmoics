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
	double theta, phi, x=1.0, y=2.0, z=3.0,r;
	double *Y;
	int Lmax, memsize_Y;
    	struct timespec t_start, t_end, t1, t2;
    	double time_secs;
    	r = sqrt(x*x+y*y+z*z);

    	my_gettime(&t1);
 
	Lmax=3;
	theta=atan(sqrt(x*x+y*y)/z);
	phi = atan(y/x);
	memsize_Y = (Lmax+1)*(Lmax+1);
	
    	my_gettime(&t_start);
    	Y = (double *) malloc(sizeof(double)*memsize_Y);
	//dY_phi = (double *) malloc(sizeof(double)*memsize_Y);
	//dY_th = (double *) malloc(sizeof(double)*memsize_Y);
	my_gettime(&t_end);
	time_secs = elapsed_time(&t_start, &t_end);
	printf("Run-time of the memory allocation: %.3lf ms\n", time_secs*1000.0);
 
	my_gettime(&t_start);
	sph_harmonics_real(theta, phi, Lmax, Y);
	my_gettime(&t_end);
	// time in seconds
	time_secs = elapsed_time(&t_start, &t_end);
	printf("Run-time of the sph_harmonics: %.3lf ms\n", time_secs*1000.0);

 
	free(Y);
	//free(dY_phi);
	//free(dY_th);

    	my_gettime(&t2);
    	time_secs = elapsed_time(&t1, &t2);
    	printf("Run-time of total program: %.3lf ms\n", time_secs*1000.0);
    	
    	double Ylm[1];
    	
    	my_gettime(&t1);
    	for (int i=0; i < Lmax; i++){
    		for (int j=-i; j<=i; j++){
    			RealSphericalHarmonic(1, &x, &y, &z, &r, i, j, Ylm);
    		}
    	}
    	my_gettime(&t2);
    	time_secs = elapsed_time(&t1, &t2);
    	printf("Run-time of total program SPARC: %.3lf ms\n", time_secs*1000.0);
    	
    	RealSphericalHarmonic(1, &x, &y, &z, &r, Lmax, Lmax, Ylm);
    	printf("SPARC: %f New: %f\n",Ylm[0], Y[memsize_Y-1]);

	return 0;
}
