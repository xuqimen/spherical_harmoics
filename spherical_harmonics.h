/**
 * @file spherical_harmonics.h
 * @author Shashikant Kumar (shashikant@gatech.edu)
 *         Qimen Xu (qimenxu@foxmail.com)
 * @brief Spherical harmonic functions.
 * @version 0.1
 * @date 2022-04-19
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef _SPHERICAL_HARMONICS_H
#define _SPHERICAL_HARMONICS_H


/*
computeP function calculates the normalized associated Legendre's polynomial. The normalization is done such that
the corresponding spherical harmonics functions are orthonormal. 

[Input]
1. L: maximum 'l' index (orbital angular momentum number)
2. A,B: coefficients for the recursion precomputed
3. x: cos(theta)
[Output]
1. P: pointer to the array containing the Legendre polynomial. Values for only m>0 are calculated.
*/
void computeP( const size_t L ,
			const double * const A , const double * const B ,
			double * const P , const double x );



/*
computeY function calculates the spherical harmonics. 

[Input]
1. L: maximum 'l' index (orbital angular momentum number)
2. P: pointer to the array containing the Legendre polynomial. Values for only m>0 are stored.
3. phi
[Output]
1. Y: pointer to the array containing the Spherical harmonics for all combinations 0 <= l <= L and -l <=m <= l. 
*/
void computeY( const size_t L , const double * const P ,
	double complex * const Y ,  const double phi );

void computeP_real( const size_t L ,
	const double * const A , const double * const B ,
	double * const P , const double x );

void computeY_real( const size_t L , const double * const P ,
	double * const Y , const double phi );

void sph_harmonics_real(const double theta, const double phi, const int LL,
				double * const Y, double * const dY_theta,
				double * const dY_phi );

/*
sph_harmonics function calculates the spherical harmonics and its derivatives with respect to theta and phi. 

[Input]
1. theta, phi: spherical coordinate of the point on sphere [ theta \in [0 \pi], phi \in [0 2\pi] ]
2. LL: maximum 'l' index (orbital angular momentum number)
[Output]
1. Y: pointer to the array containing the Spherical harmonics for all combinations 0 <= l <= L and -l <=m <= l. 
2. dY_theta: pointer to the array containing the derrivative of  Spherical harmonics w.r.t theta for all combinations 0 <= l <= L and -l <=m <= l. 
3. dY_phi: pointer to the array containing the derrivative of  Spherical harmonics w.r.t phi for all combinations 0 <= l <= L and -l <=m <= l. 
*/
void sph_harmonics(const double theta, const double phi, const int LL,
				double complex * const Y, double complex * const dY_theta,
				double complex * const dY_phi );

#endif
