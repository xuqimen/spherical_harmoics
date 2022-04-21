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
#include <assert.h>

#include "tools.h"
#include "spherical_harmonics.h"

#define PT(l1, m1) ((m1) +((l1) *((l1) +1))/2)
#define YR(l2 , m2 ) (( m2 ) +( l2 ) +(( l2 ) *( l2) ) )
#define max(x,y) ((x) > (y) ? (x) : (y))


int test() {
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


/**
 * @brief Test real spherical harmonics against explicit formula.
 * 
 * @param Lmax Maximum l (here we assume all l <= Lmax will be calculated).
 * @param n Number of coordinates to evaluate.
 */
void test_results_real(int Lmax, int n)
{
	struct timespec t_start, t_end, t1, t2;
	double time_secs;
	my_gettime(&t_start);

	// set up test input data
	// int n = 1000;
	double *x = (double *) malloc(n * sizeof(*x));
	double *y = (double *) malloc(n * sizeof(*y));
	double *z = (double *) malloc(n * sizeof(*z));
	double *r = (double *) malloc(n * sizeof(*r));
	double *theta = (double *) malloc (n * sizeof(*theta));
	double *phi = (double *) malloc (n * sizeof(*phi));
	assert(x != NULL && y != NULL && z != NULL);
	assert(r != NULL && theta != NULL && phi != NULL);

	double rand_min = -1.0, rand_max = 1.0;
	srand(1); // set up seed for random number generator
	// generate n random coordinates (x,y,z)
	SetRandMat(x, n, 1, rand_min, rand_max);
	SetRandMat(y, n, 1, rand_min, rand_max);
	SetRandMat(z, n, 1, rand_min, rand_max);

	// convert to spherical coordiantes, i.e., find radius, theta, phi
	// for (int i = 0; i < n; i++) {
	// 	r[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
	// 	theta[i] = atan(sqrt(x[i]*x[i]+y[i]*y[i])/z[i]);
	// 	phi[i] = atan(y[i]/x[i]);
	// }
	cart2polar(n, x, y, z, r, theta, phi);

	int memsize_Y = (Lmax+1)*(Lmax+1);
	double *Y = (double *) malloc(memsize_Y * n * sizeof(*Y));
	assert(Y != NULL);
	double *Ylm = (double *) malloc(n * sizeof(*Y));
	assert(Ylm != NULL);

	// use new general routine to find Y_lm for all l <= Lmax
	my_gettime(&t1); // start timer
	for (int i = 0; i < n; i++) {
		sph_harmonics_real(theta[i], phi[i], Lmax, &Y[i*memsize_Y]);
	}
	my_gettime(&t2); // stop timer
	double t_sph_harmonics = elapsed_time(&t1, &t2);

	// use explicit formula to find Ylm for all l <= Lmax
	double t_sparc = 0.0;
	double error = 0.0;
	double tol = 1e-12;
	for (int l = 0; l <= Lmax; l++) {
		for (int m = -l; m <= l; m++) {
			my_gettime(&t1); // start timer
			// use SPARC routine to find spherical harmonics for l,m
			RealSphericalHarmonic(n, x, y, z, r, l, m, Ylm);
			my_gettime(&t2); // stop timer
			t_sparc += elapsed_time(&t1, &t2);
			
			// compare results with general routine
			double error_lm = 0.0;
			error_lm = check_double_arrays(Ylm, 1, &Y[YR(l, m)], memsize_Y, n, tol);
			error = max(error, error_lm);
			printf(BLU "l = %2d, m = %2d, error_lm = %.3e\n" RESET, l, m, error_lm);
			if (error_lm > tol) {
				for (int i = 0; i < n; i++) {
					printf("(x,y,z) = (%.3f,%.3f,%.3f), l = %2d, m = %2d, Ylm = %.3f, Y[l,m] = %.3f, error_lm = %.3e\n",
							x[i],y[i],z[i],l,m,Ylm[i],Y[YR(l, m)+i*memsize_Y],error_lm);
				}
			}
		}
	}
	if (error > tol) {
		printf(RED "One or more test failed!\n" RESET);
	} else {
		printf(GRN "Success! All tests passed!\n" RESET);
	}

	free(x); free(y); free(z);
	free(r); free(theta); free(phi);
	free(Y); free(Ylm);

	my_gettime(&t_end);
	// time in seconds
	time_secs = elapsed_time(&t_start, &t_end);
	
	printf("===============\n");
	printf("= Timing info =\n");
	printf("===============\n");
	printf("Run-time of sph_harmonics: %.3lf ms\n", t_sph_harmonics*1000.0);
	printf("Run-time of SPARC routine: %.3lf ms\n", t_sparc*1000.0);
	printf("Total run-time of the test: %.3lf ms\n", time_secs*1000.0);
}

void test_results_real_vector(int Lmax, int n)
{
	struct timespec t_start, t_end, t1, t2;
	double time_secs;
	my_gettime(&t_start);

	// set up test input data
	// int n = 1000;
	double *x = (double *) malloc(n * sizeof(*x));
	double *y = (double *) malloc(n * sizeof(*y));
	double *z = (double *) malloc(n * sizeof(*z));
	double *r = (double *) malloc(n * sizeof(*r));
	double *theta = (double *) malloc (n * sizeof(*theta));
	double *phi = (double *) malloc (n * sizeof(*phi));
	assert(x != NULL && y != NULL && z != NULL);
	assert(r != NULL && theta != NULL && phi != NULL);

	double rand_min = -1.0, rand_max = 1.0;
	srand(1); // set up seed for random number generator
	// generate n random coordinates (x,y,z)
	SetRandMat(x, n, 1, rand_min, rand_max);
	SetRandMat(y, n, 1, rand_min, rand_max);
	SetRandMat(z, n, 1, rand_min, rand_max);

	// convert to spherical coordiantes, i.e., find radius, theta, phi
	// for (int i = 0; i < n; i++) {
	// 	r[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
	// 	theta[i] = atan(sqrt(x[i]*x[i]+y[i]*y[i])/z[i]);
	// 	phi[i] = atan(y[i]/x[i]);
	// }
	cart2polar(n, x, y, z, r, theta, phi);

	int memsize_Y = (Lmax+1)*(Lmax+1);
	double *Y = (double *) malloc(memsize_Y * n * sizeof(*Y));
	assert(Y != NULL);
	double *Ylm = (double *) malloc(n * sizeof(*Y));
	assert(Ylm != NULL);

	// use new general routine to find Y_lm for all l <= Lmax
	my_gettime(&t1); // start timer
	sph_harmonics_real_vector(n, theta, phi, Lmax, Y);
	my_gettime(&t2); // stop timer
	double t_sph_harmonics = elapsed_time(&t1, &t2);

	// use explicit formula to find Ylm for all l <= Lmax
	double t_sparc = 0.0;
	double error = 0.0;
	double tol = 1e-12;
	for (int l = 0; l <= Lmax; l++) {
		for (int m = -l; m <= l; m++) {
			my_gettime(&t1); // start timer
			// use SPARC routine to find spherical harmonics for l,m
			RealSphericalHarmonic(n, x, y, z, r, l, m, Ylm);
			my_gettime(&t2); // stop timer
			t_sparc += elapsed_time(&t1, &t2);
			
			// compare results with general routine
			double error_lm = 0.0;
			error_lm = check_double_arrays(Ylm, 1, &Y[YR(l, m)], memsize_Y, n, tol);
			error = max(error, error_lm);
			printf(BLU "l = %2d, m = %2d, error_lm = %.3e\n" RESET, l, m, error_lm);
			if (error_lm > tol) {
				for (int i = 0; i < n; i++) {
					printf("(x,y,z) = (%.3f,%.3f,%.3f), l = %2d, m = %2d, Ylm = %.3f, Y[l,m] = %.3f, error_lm = %.3e\n",
							x[i],y[i],z[i],l,m,Ylm[i],Y[YR(l, m)+i*memsize_Y],error_lm);
				}
			}
		}
	}
	if (error > tol) {
		printf(RED "One or more test failed!\n" RESET);
	} else {
		printf(GRN "Success! All tests passed!\n" RESET);
	}

	free(x); free(y); free(z);
	free(r); free(theta); free(phi);
	free(Y); free(Ylm);

	my_gettime(&t_end);
	// time in seconds
	time_secs = elapsed_time(&t_start, &t_end);
	
	printf("===============\n");
	printf("= Timing info =\n");
	printf("===============\n");
	printf("Run-time of sph_harmonics: %.3lf ms\n", t_sph_harmonics*1000.0);
	printf("Run-time of SPARC routine: %.3lf ms\n", t_sparc*1000.0);
	printf("Total run-time of the test: %.3lf ms\n", time_secs*1000.0);
}



int main(int argc, char *argv[]) {
	// test();

	int Lmax = 6, n = 100;
	if (argc == 3) {
		Lmax = atoi(argv[1]);
		n = atoi(argv[2]);
	}
	test_results_real_vector(Lmax, n);

	return 0;
}
