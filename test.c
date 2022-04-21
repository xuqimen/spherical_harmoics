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
#define YR(l2 , m2 ) (( m2 ) +( l2 ) +(( l2 ) *( l2) ) )
int main() {
	double theta, phi, x=-1.0, y=-1.0, z=-1.0, r;
	double *Y;
	int Lmax, memsize_Y;
    struct timespec t_start, t_end, t1, t2;
    double time_secs;
    r = sqrt(x*x+y*y+z*z);

    my_gettime(&t1);
 
	Lmax=3;
	const double PI=3.141592653589793;
	theta=atan(sqrt(x*x+y*y)/z);
	if (theta<0) theta = PI + theta;

	phi = atan(y/x);
	if (x<0) phi = PI+phi;  // 2nd Quandrant
	if (x>0 && y<0) phi = 2*PI+phi;  // 4th Quadrant
	if (x==0 && y <0) phi = 3*PI/2;
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

	
	for (int i=0; i < Lmax; i++){
		for (int j=-i; j<=i; j++){
			RealSphericalHarmonic(1, &x, &y, &z, &r, i, j, Ylm);
			printf("(x,y,z)=(%f,%f,%f), l=%d, m=%d, SPARC: %f New: %f, error: %10.9f\n",x,y,z,i,j, Ylm[0], Y[YR(i,j)], Ylm[0]-Y[YR(i,j)]);
		}
	}
	


	free(Y);
	return 0;
}
