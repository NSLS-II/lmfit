/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     hat.c
 *
 * Contents: An example with m=2, n=1 from Madsen et al. 2004, Example 3.3
 *
 *           Minimize the norm of a 2-dimensional vectorial function
 *           of a 1-dimensional parameter: F(p) = ( p+1, lambda*p^2+p-1 ).
 *           So, ||F|| is kind of an asymetric mexican hat function.
 *
 *           There is a stationary point at par = 0. Madsen et al. indicate
 *           incorrectly that it is a global minimum for lambda<1.
 *           It is so only for lambda below about 0.5. For higher values
 *           of lambda, this example allows to study the dependence of
 *           fit results on starting values.
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  The Mexican-Hat-Do-What-You-Want-License Rev. 7.31
 */
 
#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>


void evaluate_m2n1( const double *par, int m_dat, const void *data,
                    double *fvec, int *info )
{
    fvec[0] = par[0]+1;
    fvec[1] = (*((double*)data))*par[0]*par[0] + par[0] - 1;
}


int main( int argc, char **argv )
{
    /* parameter lambda */

    if( argc!=3 ){
        fprintf( stderr, "usage: m2n1 lambda p_start\n" );
        exit(-1);
    }
    double lambda = atof( argv[1] );

    /* parameter vector */

    int n_par = 1; // number of parameters in model function f
    double par[1] = { atof( argv[2] ) }; // arbitrary starting value

    /* data points */

    int m_dat = 2;

    /* auxiliary parameters */

    lm_status_struct status; // to receive status information
    lm_control_struct control = lm_control_double;
    control.printflags = 3; // monitor status (+1) and parameters (+2)

    /* perform the fit */

    printf( "Fitting:\n" );
    lmmin( n_par, par, m_dat, (const void*) &lambda,
           evaluate_m2n1, &control, &status, lm_printout_std );

    /* print results */

    printf( "\nResults:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.info] );

    printf("obtained parameters:\n");
    int i;
    for ( i=0; i<n_par; ++i )
	printf("  par[%i] = %19.11f\n", i, par[i]);
    printf("obtained norm:\n  %19.11f\n", status.fnorm );

    return 0;
}
