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
	P[PT(0,0)] = temp ;
	if ( L > 0) {
		const double SQRT3 = 1.732050807568877;
		P[PT(1,0)] = x * SQRT3 * temp ;
		const double SQRT3DIV2 = -1.2247448713915890491;
		temp = SQRT3DIV2 * sintheta * temp ;
		P[PT(1,1)] = temp ;

		for ( size_t l =2; l <= L ; l ++) {
			for ( size_t m =0; m <l -1; m ++) {
				P[PT(l,m)] = A[PT(l,m)]*( x * P[PT(l-1,m)]
				+ B[PT(l,m)]* P[PT(l-2,m)]) ;
			}
			P[PT(l,l-1)] = x * sqrt(2*( l-1) +3) * temp ;
			temp = - sqrt(1.0+0.5/l) * sintheta * temp ;
			P[PT(l,l)] = temp ;
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
	double c1 = 1.0 , c2 = cos(phi) ;
	double s1 = 0.0 , s2 = - sin(phi) ;
	double tc = 2.0 * c2 ;
	for ( size_t m =1; m <= L ; m ++) {
		double s = tc * s1 - s2 ;
		double c = tc * c1 - c2 ;
		s2 = s1 ; s1 = s ; c2 = c1 ; c1 = c ;
		for ( size_t l = m ; l <= L ; l ++) {
			
			temp1 = P[PT(l,m)] * c + P[ PT(l,m) ] * s * I ;
			if (m%2 == 0){
				temp2 =  P[PT(l,m)] * c - P[PT(l,m)] * s * I ;
			} else {
				temp2 =  -P[PT(l,m)] * c + P[PT(l,m)] * s * I ;
			}
			
			Y[YR(l, m)] = temp1 ;
			Y[YR(l, -m)] = temp2 ;

		}
	}
}

void computeP_real( const size_t L ,
	const double * const A , const double * const B ,
	double * const P , const double x )
{
	const double sintheta = sqrt (1.0 - x*x) ;
	double temp = 0.282094791773878 ; // = sqrt (1/ 4 M_PI )
	P[PT(0, 0) ] = temp ;
	if (L > 0) {
		const double SQRT3 = 1.7320508075688772935 ;
		P[PT(1, 0) ] = x * SQRT3 * temp ;
		const double SQRT3DIV2 = -1.2247448713915890491;
		temp = SQRT3DIV2 * sintheta * temp ;
		P[PT(1, 1)] = temp ;

		for ( size_t l = 2; l <= L ; l ++) {
			for ( size_t m = 0; m < l-1; m ++) {
				P[PT(l,m)] = A[PT(l, m)]*(x*P[PT(l-1, m )]
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
		Y[YR(l, 0)] = P[PT(l, 0)] ;
	
	const double SQRT2 = 1.414213562373095 ;
	
	double c1 = 1.0 , c2 = cos(phi);
	double s1 = 0.0 , s2 = -sin(phi) ;
	double tc = 2.0 * c2 ;
	for (size_t m = 1; m <= L ; m ++) {
		double s = tc * s1 - s2;
		double c = tc * c1 - c2;
		s2 = s1 ; s1 = s ; c2 = c1 ; c1 = c ;
		for ( size_t l = m ; l <= L ; l ++) {
			if (m%2==0){
				Y[YR(l, - m)] = P[PT(l, m)] * s*SQRT2 ;
				Y[YR(l, m)] = P[PT(l, m)] * c*SQRT2;
			} else {
				Y[YR(l, - m)] = -P[PT(l, m)] * s*SQRT2 ;
				Y[YR(l, m)] = -P[PT(l, m)] * c*SQRT2;
			}
			
		}
	}
}

void sph_harmonics_real(const double theta, const double phi, const int LL,
				double * const Y) {
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
			A[PT(l, m)] = sqrt((4* ls -1.0) /( ls - ms ) ) ;
			B[PT(l, m)] = - sqrt(( lm1s - ms ) /(4* lm1s -1.0) ) ;
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
			A[PT(l, m)] = sqrt((4* ls -1.0) /( ls - ms ) ) ;
			B[PT(l, m)] = - sqrt(( lm1s - ms ) /(4* lm1s -1.0) ) ;
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

void RealSphericalHarmonic(const int len, double *x, double *y,double *z, double *r, const int l, const int m, double *Ylm)
{
    // only l=0,1,2,3,4,5,6 implemented for now

    //double pi=M_PI;
    double p;                   	  
    int i; 
    
    /* l = 0 */
    double C00 = 0.282094791773878; // 0.5*sqrt(1/pi)
    /* l = 1 */
    double C1m1 = 0.488602511902920; // sqrt(3/(4*pi))
    double C10 = 0.488602511902920; // sqrt(3/(4*pi))
    double C1p1 = 0.488602511902920; // sqrt(3/(4*pi))
    /* l = 2 */
    double C2m2 = 1.092548430592079; // 0.5*sqrt(15/pi)
    double C2m1 = 1.092548430592079; // 0.5*sqrt(15/pi)  
    double C20 =  0.315391565252520; // 0.25*sqrt(5/pi)
    double C2p1 = 1.092548430592079; // 0.5*sqrt(15/pi)  
    double C2p2 = 0.546274215296040; // 0.25*sqrt(15/pi)
    /* l = 3 */
    double C3m3 = 0.590043589926644; // 0.25*sqrt(35/(2*pi))   
    double C3m2 = 2.890611442640554; // 0.5*sqrt(105/(pi))
    double C3m1 = 0.457045799464466; // 0.25*sqrt(21/(2*pi))
    double C30 =  0.373176332590115; // 0.25*sqrt(7/pi)
    double C3p1 = 0.457045799464466; // 0.25*sqrt(21/(2*pi))
    double C3p2 = 1.445305721320277; //  0.25*sqrt(105/(pi))
    double C3p3 = 0.590043589926644; //  0.25*sqrt(35/(2*pi))
    /* l = 4 */
    double C4m4 = 2.503342941796705; // (3.0/4.0)*sqrt(35.0/pi)
    double C4m3 = 1.770130769779930; // (3.0/4.0)*sqrt(35.0/(2.0*pi))   
    double C4m2 = 0.946174695757560; // (3.0/4.0)*sqrt(5.0/pi)
    double C4m1 = 0.669046543557289; // (3.0/4.0)*sqrt(5.0/(2.0*pi))
    double C40 =  0.105785546915204; // (3.0/16.0)*sqrt(1.0/pi)
    double C4p1 = 0.669046543557289; // (3.0/4.0)*sqrt(5.0/(2.0*pi))
    double C4p2 = 0.473087347878780; // (3.0/8.0)*sqrt(5.0/(pi))
    double C4p3 = 1.770130769779930; // (3.0/4.0)*sqrt(35.0/(2.0*pi))
    double C4p4 = 0.625835735449176; // (3.0/16.0)*sqrt(35.0/pi) 
    /* l = 5 */
    double C5m5 = 0.656382056840170; // (3.0*sqrt(2.0*77.0/pi)/32.0)
    double C5m4 = 2.075662314881042; // (3.0/16.0)*sqrt(385.0/pi)
    double C5m3 = 0.489238299435250; // (sqrt(2.0*385.0/pi)/32.0)
    double C5m2 = 4.793536784973324; // (1.0/8.0)*sqrt(1155.0/pi)*2.0
    double C5m1 = 0.452946651195697; // (1.0/16.0)*sqrt(165.0/pi)
    double C50 =  0.116950322453424; // (1.0/16.0)*sqrt(11.0/pi)
    double C5p1 = 0.452946651195697; // (1.0/16.0)*sqrt(165.0/pi)
    double C5p2 = 2.396768392486662; // (1.0/8.0)*sqrt(1155.0/pi)
    double C5p3 = 0.489238299435250; // (sqrt(2.0*385.0/pi)/32.0)
    double C5p4 = 2.075662314881042; // (3.0/16.0)*sqrt(385.0/pi)
    double C5p5 = 0.656382056840170; // (3.0*sqrt(2.0)/32.0)*sqrt(77.0/pi)
    /* l = 6 */
    double C6m6 = 0.683184105191914; // (sqrt(2.0*3003.0/pi)/64.0)
    double C6m5 = 2.366619162231752; // (3.0/32.0)*sqrt(2.0*1001.0/pi)
    double C6m4 = 0.504564900728724; // (3.0/32.0)*sqrt(91.0/pi)
    double C6m3 = 0.921205259514923; // (sqrt(2.0*1365.0/pi)/32.0)
    double C6m2 = 0.460602629757462; // (sqrt(2.0*1365/pi)/64.0)
    double C6m1 = 0.582621362518731; // (sqrt(273.0/pi)/16.0)
    double C60 =  0.0635692022676284;// (sqrt(13.0/pi)/32.0)
    double C6p1 = 0.582621362518731; // (sqrt(273.0/pi)/16.0)
    double C6p2 = 0.460602629757462; // (sqrt(2.0*1365.0/pi)/64.0)
    double C6p3 = 0.921205259514923; // (sqrt(2.0*1365.0/pi)/32.0)
    double C6p4 = 0.504564900728724; // (3.0/32.0)*sqrt(91.0/pi)
    double C6p5 = 2.366619162231752; // (3.0/32.0)*sqrt(2.0*1001.0/pi)
    double C6p6 = 0.683184105191914; // (sqrt(2.0*3003.0/pi)/64.0)
    
    switch (l)
    {
        /* l = 0 */
        case 0:
            for (i = 0; i < len; i++) Ylm[i] = C00;
            break;
        
        /* l = 1 */
        case 1: 
            switch (m) 
            {
                case -1: /* m = -1 */
                    for (i = 0; i < len; i++) 
                        Ylm[i] = C1m1 * (y[i] / r[i]);
                    break;
                
                case 0: /* m = 0 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C10 * (z[i] / r[i]);
                    break;
                
                case 1: /* m = 1 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C1p1 * (x[i] / r[i]);
                    break;
                
                /* incorrect m */
                default: printf("<m> must be an integer between %d and %d!\n", -l, l); break;
            }
            break;
        
        /* l = 2 */
        case 2: 
            switch (m) 
            {      
                case -2: /* m = -2 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C2m2 * (x[i]*y[i])/(r[i]*r[i]);
                    break;

                case -1: /* m = -1 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C2m1*(y[i]*z[i])/(r[i]*r[i]);
                    break;

                case 0: /* m = 0 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C20*(-x[i]*x[i] - y[i]*y[i] + 2.0*z[i]*z[i])/(r[i]*r[i]);
                    break;
                
                case 1: /* m = 1 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C2p1*(z[i]*x[i])/(r[i]*r[i]);
                    break;
                
                case 2: /* m = 2 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C2p2*(x[i]*x[i] - y[i]*y[i])/(r[i]*r[i]);
                    break;
                
                /* incorrect m */
                default: printf("<m> must be an integer between %d and %d!\n", -l, l); 
                         break;
            }
            break;

        /* l = 3 */
        case 3: 
            switch (m) 
            {   
                case -3: /* m = -3 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C3m3*(3*x[i]*x[i] - y[i]*y[i])*y[i]/(r[i]*r[i]*r[i]);
                    break;
               
                case -2: /* m = -2 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C3m2*(x[i]*y[i]*z[i])/(r[i]*r[i]*r[i]);
                    break;

                case -1: /* m = -1 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C3m1*y[i]*(4*z[i]*z[i] - x[i]*x[i] - y[i]*y[i])/(r[i]*r[i]*r[i]);
                    break;

                case 0: /* m = 0 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C30*z[i]*(2*z[i]*z[i]-3*x[i]*x[i]-3*y[i]*y[i])/(r[i]*r[i]*r[i]);
                    break;
                
                case 1: /* m = 1 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C3p1*x[i]*(4*z[i]*z[i] - x[i]*x[i] - y[i]*y[i])/(r[i]*r[i]*r[i]);
                    break;
                
                case 2: /* m = 2 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C3p2*z[i]*(x[i]*x[i] - y[i]*y[i])/(r[i]*r[i]*r[i]);
                    break;
                
                case 3: /* m = 3 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C3p3*x[i]*(x[i]*x[i]-3*y[i]*y[i])/(r[i]*r[i]*r[i]);
                    break;
                
                /* incorrect m */
                default: printf("<m> must be an integer between %d and %d!\n", -l, l); 
                         break;
            }
            break;
        
        /* l = 4 */
        case 4: 
            switch (m) 
            {     
                case -4: /* m = -4 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4m4*(x[i]*y[i]*(x[i]*x[i]-y[i]*y[i]))/(r[i]*r[i]*r[i]*r[i]);
                    break;
                
                case -3: /* m = -3 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4m3*(3.0*x[i]*x[i]-y[i]*y[i])*y[i]*z[i]/(r[i]*r[i]*r[i]*r[i]);
                    break;
               
                case -2: /* m = -2 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4m2*x[i]*y[i]*(7.0*z[i]*z[i]-r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]);
                    break;

                case -1: /* m = -1 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4m1*y[i]*z[i]*(7.0*z[i]*z[i]-3.0*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]);
                    break;

                case 0: /* m = 0 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C40*(35.0*z[i]*z[i]*z[i]*z[i]-30.0*z[i]*z[i]*r[i]*r[i]+3.0*r[i]*r[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]);
                    break;
                
                case 1: /* m = 1 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4p1*x[i]*z[i]*(7.0*z[i]*z[i]-3.0*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]);
                    break;
                
                case 2: /* m = 2 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4p2*(x[i]*x[i]-y[i]*y[i])*(7.0*z[i]*z[i]-r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]);
                    break;
                
                case 3: /* m = 3 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4p3*(x[i]*x[i]-3.0*y[i]*y[i])*x[i]*z[i]/(r[i]*r[i]*r[i]*r[i]);
                    break;
                
                case 4: /* m = 4 */
                    for (i = 0; i < len; i++)
                        Ylm[i]=C4p4*(x[i]*x[i]*(x[i]*x[i]-3.0*y[i]*y[i]) - y[i]*y[i]*(3.0*x[i]*x[i]-y[i]*y[i]))/(r[i]*r[i]*r[i]*r[i]);
                    break;
                    
                /* incorrect m */
                default: printf("<m> must be an integer between %d and %d!\n", -l, l); 
                         break;
            }
            break;
      
        /* l = 5 */
        case 5: 
            //p = sqrt(x[i]*x[i]+y[i]*y[i]);
            switch (m) 
            {   
                case -5: /* m = -5 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5m5*(8.0*x[i]*x[i]*x[i]*x[i]*y[i]-4.0*x[i]*x[i]*y[i]*y[i]*y[i] + 4.0*pow(y[i],5)-3.0*y[i]*p*p*p*p)/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case -4: /* m = -4 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5m4*(4.0*x[i]*x[i]*x[i]*y[i] - 4.0*x[i]*y[i]*y[i]*y[i])*z[i]/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case -3: /* m = -3 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5m3*(3.0*y[i]*p*p - 4.0*y[i]*y[i]*y[i])*(9.0*z[i]*z[i]-r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
               
                case -2: /* m = -2 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5m2*x[i]*y[i]*(3.0*z[i]*z[i]*z[i]-z[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;

                case -1: /* m = -1 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5m1*y[i]*(21.0*z[i]*z[i]*z[i]*z[i] - 14.0*r[i]*r[i]*z[i]*z[i]+r[i]*r[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;

                case 0: /* m = 0 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C50*(63.0*z[i]*z[i]*z[i]*z[i]*z[i] -70.0*z[i]*z[i]*z[i]*r[i]*r[i] + 15.0*z[i]*r[i]*r[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 1: /* m = 1 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5p1*x[i]*(21.0*z[i]*z[i]*z[i]*z[i] - 14.0*r[i]*r[i]*z[i]*z[i]+r[i]*r[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 2: /* m = 2 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5p2*(x[i]*x[i]-y[i]*y[i])*(3.0*z[i]*z[i]*z[i] - r[i]*r[i]*z[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 3: /* m = 3 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5p3*(4.0*x[i]*x[i]*x[i]-3.0*p*p*x[i])*(9.0*z[i]*z[i]-r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 4: /* m = 4 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5p4*(4.0*(x[i]*x[i]*x[i]*x[i]+y[i]*y[i]*y[i]*y[i])-3.0*p*p*p*p)*z[i]/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 5: /* m = 5 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C5p5*(4.0*x[i]*x[i]*x[i]*x[i]*x[i] + 8.0*x[i]*y[i]*y[i]*y[i]*y[i] -4.0*x[i]*x[i]*x[i]*y[i]*y[i] -3.0*x[i]*p*p*p*p)/(r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                    
                /* incorrect m */
                default: printf("<m> must be an integer between %d and %d!\n", -l, l); 
                         break;
            }
            break;

        /* l = 6 */
        case 6: 
            //p = sqrt(x[i]*x[i]+y[i]*y[i]);
            switch (m) 
            {   
                case -6: /* m = -6 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6m6*(12.0*pow(x[i],5)*y[i]+12.0*x[i]*pow(y[i],5) - 8.0*x[i]*x[i]*x[i]*y[i]*y[i]*y[i]-6.0*x[i]*y[i]*pow(p,4))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case -5: /* m = -5 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6m5*(8.0*pow(x[i],4)*y[i] - 4.0*x[i]*x[i]*y[i]*y[i]*y[i] + 4.0*pow(y[i],5) -3.0*y[i]*pow(p,4))*z[i]/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case -4: /* m = -4 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C6m4*(4.0*x[i]*x[i]*x[i]*y[i] -4.0*x[i]*y[i]*y[i]*y[i])*(11.0*z[i]*z[i]-r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    break;
                
                case -3: /* m = -3 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6m3*(-4.0*y[i]*y[i]*y[i] + 3.0*y[i]*p*p)*(11.0*z[i]*z[i]*z[i] - 3.0*z[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
               
                case -2: /* m = -2 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6m2*(2.0*x[i]*y[i])*(33.0*pow(z[i],4)-18.0*z[i]*z[i]*r[i]*r[i] + pow(r[i],4))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;

                case -1: /* m = -1 */
                    for (i = 0; i < len; i++)
                        Ylm[i] = C6m1*y[i]*(33.0*pow(z[i],5)-30.0*z[i]*z[i]*z[i]*r[i]*r[i] +5.0*z[i]*pow(r[i],4))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    break;

                case 0: /* m = 0 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C60*(231.0*pow(z[i],6)-315*pow(z[i],4)*r[i]*r[i] + 105.0*z[i]*z[i]*pow(r[i],4) -5.0*pow(r[i],6))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 1: /* m = 1 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6p1*x[i]*(33.0*pow(z[i],5)-30.0*z[i]*z[i]*z[i]*r[i]*r[i] +5.0*z[i]*pow(r[i],4))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 2: /* m = 2 */
                    for (i = 0; i < len; i++) {
                        //p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6p2*(x[i]*x[i]-y[i]*y[i])*(33.0*pow(z[i],4) - 18.0*z[i]*z[i]*r[i]*r[i] + pow(r[i],4))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 3: /* m = 3 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6p3*(4.0*x[i]*x[i]*x[i] -3.0*x[i]*p*p)*(11.0*z[i]*z[i]*z[i] - 3.0*z[i]*r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 4: /* m = 4 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6p4*(4.0*pow(x[i],4)+4.0*pow(y[i],4) -3.0*pow(p,4))*(11.0*z[i]*z[i] -r[i]*r[i])/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 5: /* m = 5 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6p5*(4.0*pow(x[i],5) + 8.0*x[i]*pow(y[i],4)-4.0*x[i]*x[i]*x[i]*y[i]*y[i]-3.0*x[i]*pow(p,4))*z[i]/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                
                case 6: /* m = 6 */
                    for (i = 0; i < len; i++) {
                        p = sqrt(x[i]*x[i]+y[i]*y[i]);
                        Ylm[i] = C6p6*(4.0*pow(x[i],6)-4.0*pow(y[i],6) +12.0*x[i]*x[i]*pow(y[i],4)-12.0*pow(x[i],4)*y[i]*y[i] + 3.0*y[i]*y[i]*pow(p,4)-3.0*x[i]*x[i]*pow(p,4))/(r[i]*r[i]*r[i]*r[i]*r[i]*r[i]);
                    }
                    break;
                    
                /* incorrect m */
                default: printf("<m> must be an integer between %d and %d!\n", -l, l); 
                         break;
            }
            break;
        
        default: printf("<l> must be an integer between 0 and 6!\n"); break;
    }
    
    if (l > 0) {
        for (i = 0; i < len; i++) {
            if (r[i] < 1e-10) Ylm[i] = 0.0;
        }
    }
}

