/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     morobropro.c
 *
 * Contents: Modified Rosenbrock Problem,
 *           according to Madsen et al. 2004, Example 3.13
 *
 *           Minimize the norm of a 3-dimensional vectorial function
 *           of a 2-dimensional parameter vector:
 *           F(p) = ( 10(p1-p0^2), 1-p0, lambda ).
 *
 *           The analytical solution is par = ( 1, 1 ),
 *           independent of the external parameter lambda.
 *
 *           Madsen et al. report that conventional LM fails for
 *           lambda=1e4. Fortunately, this is not the case here.
 *           On a x86, perfect results are found up to lambda>2e6.
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 
#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>


void evaluate_morobropro( const double *par, int m_dat, const void *data,
                          double *fvec, int *info )
{
    fvec[0] = 10*(par[1]-par[0]*par[0]);
    fvec[1] = 1 - par[0];
    fvec[2] = *((double*)data);
}


int main( int argc, char **argv )
{
    /* parameter lambda */

    if( argc!=2 ){
        fprintf( stderr, "usage: morobropro lambda\n" );
        exit(-1);
    }
    double lambda = atof( argv[1] );

    /* parameter vector */

    int n_par = 2; // number of parameters in model function f
    double par[2] = { -1.2, 1 }; // arbitrary starting value

    /* data points */

    int m_dat = 3;

    /* auxiliary parameters */

    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    control.printflags = 3; // monitor status (+1) and parameters (+2)

    /* perform the fit */

    printf( "Fitting:\n" );
    lmmin( n_par, par, m_dat, (const void*) &lambda,
           evaluate_morobropro, &control, &status, lm_printout_std );

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
