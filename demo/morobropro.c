/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     morobropro.c
 *
 * Contents: Modified Rosenbrock Problem,
 *           according to Madsen et al. 2004, Example 3.13
 *
 * Author:   Joachim Wuttke 2010, following a suggestion by Mario Rudolphi
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Creative Commons Attribution Share Alike.
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
    if( argc!=1 ){
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

    lm_status_struct status; // to receive status information
    int printflags = 3;      // monitor status (+1) and parameters (+2)

    /* perform the fit */

    printf( "Fitting:\n" );
    lmmin( n_par, par, m_dat, (const void*) &lambda,
           evaluate_surface, &lm_limits_double, &status,
           lm_printout_std, printflags );

    /* print results */

    printf( "\nResults:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.info] );

    printf("obtained parameters:\n");
    int i;
    for ( i=0; i<n_par; ++i )
	printf("  par[%i] = %12g\n", i, par[i]);
    printf("obtained norm:\n  %12g\n", status.fnorm );

    return 0;
}
