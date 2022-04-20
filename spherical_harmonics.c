/* 
Implementation of Spherical Harmonics and its derivatives with respect to Theta
and Phi as implemented in the book "Numerical recipes, 3rd ed" and in the following paper:
https://arxiv.org/abs/1410.1748v1
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "tools.h"
#include "spherical_harmonics.h"

#define PT(l1, m1) ((m1) +((l1) *((l1) +1))/2)
#define YR(l2 , m2 ) (( m2 ) +( l2 ) +(( l2 ) *( l2) ) )

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
			double * const P , const double x ) {

	const double sintheta = sqrt (1.0 - x * x ) ;
	double temp = 0.282094791773878 ; // = sqrt (1/ 4 M_PI )
	P[ PT(0,0) ] = temp ;
	if ( L > 0) {
		const double SQRT3 = 1.732050807568877;
		P[ PT(1,0) ] = x * SQRT3 * temp ;
		const double SQRT3DIV2 = -1.2247448713915890491;
		temp = SQRT3DIV2 * sintheta * temp ;
		P[ PT(1,1) ] = temp ;

		for ( size_t l =2; l <= L ; l ++) {
			for ( size_t m =0; m <l -1; m ++) {
				P[ PT(l,m) ] = A[ PT(l,m) ]*( x * P[ PT(l-1,m) ]
				+ B[ PT(l,m) ]* P[ PT(l-2,m) ]) ;
			}
			P[ PT(l,l-1) ] = x * sqrt (2*( l-1) +3) * temp ;
			temp = - sqrt(1.0+0.5/l) * sintheta * temp ;
			P[ PT(l,l) ] = temp ;
		}
	}
}



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
	double complex * const Y ,  const double phi ) {

	for ( size_t l =0; l <= L ; l ++){
		Y[ YR(l,0) ] = P[ PT(l, 0) ]  ;
	}
	
	double complex temp1, temp2;
	double c1 = 1.0 , c2 = cos ( phi ) ;
	double s1 = 0.0 , s2 = - sin ( phi ) ;
	double tc = 2.0 * c2 ;
	for ( size_t m =1; m <= L ; m ++) {
		double s = tc * s1 - s2 ;
		double c = tc * c1 - c2 ;
		s2 = s1 ; s1 = s ; c2 = c1 ; c1 = c ;
		for ( size_t l = m ; l <= L ; l ++) {
			
			temp1 = P[ PT(l,m) ] * c + P[ PT(l,m) ] * s * I ;
			if (m%2 == 0){
				temp2 =  P[ PT(l,m) ] * c - P[ PT(l,m) ] * s * I ;
			} else {
				temp2 =  -P[ PT(l,m) ] * c + P[ PT(l,m) ] * s * I ;
			}
			
			Y[ YR(l, m) ] = temp1 ;
			Y[ YR(l, -m) ] = temp2 ;

		}
	}
}

void computeP_real( const size_t L ,
	const double * const A , const double * const B ,
	double * const P , const double x )
{
	const double sintheta = sqrt (1.0 - x*x) ;
	double temp = 0.39894228040143267794 ; // = sqrt (0.5/ M_PI )
	P[PT(0, 0) ] = temp ;
	if (L > 0) {
		const double SQRT3 = 1.7320508075688772935 ;
		P [PT(1, 0) ] = x * SQRT3 * temp ;
		const double SQRT3DIV2 = -1.2247448713915890491;
		temp = SQRT3DIV2 * sintheta * temp ;
		P [PT(1, 1)] = temp ;

		for ( size_t l = 2; l <= L ; l ++) {
			for ( size_t m = 0; m < l-1; m ++) {
				P [PT(l,m)] = A[PT(l, m)]*(x*P[PT(l-1, m )]
				+ B[PT(l, m )]* P[PT(l-2, m)]) ;
			}
			P[PT(l, l-1)] = x*sqrt(2*(l-1)+3)*temp ;
			temp = - sqrt(1.0+0.5/l) * sintheta * temp ;
			P[PT(l, l)] = temp ;
		}
	}
}


void computeY_real( const size_t L , const double * const P ,
	double * const Y , const double phi ) {
	for (size_t l = 0; l <= L ; l ++)
		Y[YR(l, 0)] = P[PT(l, 0)] * 0.5 * M_SQRT2 ;

	double c1 = 1.0 , c2 = cos(phi);
	double s1 = 0.0 , s2 = -sin(phi) ;
	double tc = 2.0 * c2 ;
	for (size_t m = 1; m <= L ; m ++) {
		double s = tc * s1 - s2;
		double c = tc * c1 - c2;
		s2 = s1 ; s1 = s ; c2 = c1 ; c1 = c ;
		for ( size_t l = m ; l <= L ; l ++) {
			Y[YR(l, - m)] = P[PT(l, m)] * s ;
			Y[YR(l, m)] = P[PT(l, m)] * c ;
		}
	}
}

void sph_harmonics_real(const double theta, const double phi, const int LL,
				double * const Y, double * const dY_theta,
				double * const dY_phi ) {
	int memsize_A_B_P,  memsize_Y;
	memsize_A_B_P = ((LL+1)*(LL+2))/2;

	double A[memsize_A_B_P], B[memsize_A_B_P];
	double *P;
	double x = cos(theta);
	
	memsize_Y = (LL+1)*(LL+1) ;

	P = (double *) malloc(sizeof(double)*memsize_A_B_P);
	for ( size_t l =2; l <= LL ; l++) {
		double ls = l *l , lm1s = (l -1) *( l -1) ;
		for ( size_t m =0; m <l -1; m ++) {
			double ms = m * m ;
			A[ PT(l, m) ] = sqrt ((4* ls -1.0) /( ls - ms ) ) ;
			B[ PT(l, m) ] = - sqrt (( lm1s - ms ) /(4* lm1s -1.0) ) ;
		}
	}
	
	
	computeP_real( LL , &A[0] , &B[0] , P , x );
	computeY_real( LL , P , Y ,  phi );
	
	/* Derivatives w.r.t theta and phi remaining */

	free(P);

}


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
				double complex * const dY_phi ) {
  	
	// int memsize_A_B_P,  memsize_Y;
	int memsize_A_B_P;
	memsize_A_B_P = ((LL+1)*(LL+2))/2;

	double A[memsize_A_B_P], B[memsize_A_B_P];
	double *P;
	double x = cos(theta);
	
	// memsize_Y = (LL+1)*(LL+1) ;

	P = (double *) malloc(sizeof(double)*memsize_A_B_P);
	for ( size_t l =2; l <= LL ; l++) {
		double ls = l *l , lm1s = (l -1) *( l -1) ;
		for ( size_t m =0; m <l -1; m ++) {
			double ms = m * m ;
			A[ PT(l, m) ] = sqrt ((4* ls -1.0) /( ls - ms ) ) ;
			B[ PT(l, m) ] = - sqrt (( lm1s - ms ) /(4* lm1s -1.0) ) ;
		}

	}
	
	
	computeP( LL , &A[0] , &B[0] , P , x );
	computeY( LL , P , Y ,  phi );
	
	double constant;
	constant=1;
	if (theta<0.0) constant = -1;
	double theta_abs;
	theta_abs = fabs(theta);
	dY_phi[ YR(0, 0) ] = 0.0 ;
	dY_theta[ YR(0, 0) ] = 0.0;

	double theta_cut = 0.0000000001;
	double c = cos(theta_abs), s = sin(theta_abs), ct=0.0;
	double complex emph = cos(phi) - sin(phi) *I;
	if (fabs(theta_abs) > theta_cut){
		ct = c/s;
	}
	
	for ( int l =1; l <= LL ; l++) {
		for ( int m =-l; m <=l; m++) {
			dY_phi[ YR(l, m) ] = Y[YR(l, m)] * m * I;
			if (fabs(theta_abs) > theta_cut){
				if (m<l){
					dY_theta[ YR(l, m) ] = constant*(m*ct*Y[ YR(l,m) ] + sqrt((l-m)*(l+m+1)) * emph * Y[ YR(l,m+1)]);		
				} else {
					dY_theta[ YR(l, m) ] = constant*(m*ct*Y[ YR(l,m) ]);
				}
			} else {
				dY_theta[ YR(l, m) ] = 0.0;
			}
			
		}
	} 
	

	free(P);

}


