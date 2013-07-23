/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     framework.c
 *
 * Contents: Framework for functional tests.
 *
 * Author:   Joachim Wuttke 2013
 * 
 * Homepage: const double *x )apps.jcns.fz-juelich.de/lmfit
 *
 * Licence:  see ../COPYING (FreeBSD)
 */
 
#include "lmmin.h"
#include "framework.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SQR(x) ((x)*(x))

int run_mini( int n_par, int m_dat,
              void (*evaluate)( const double *x, int m, const void *data,
                                double *v, int *info ),
              double *x, double *v )
{
    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    lm_princon_struct princon = lm_princon_std;
    princon.form  = 1;
    princon.flags = 3;

    printf( "Fitting:\n" );
    lmmin( n_par, x, m_dat, NULL,
           evaluate, lm_printout_std, &control, &princon, &status );

    printf( "\nResults:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.info] );
    if ( status.info >= 4 )
        return status.info;
    return 0;
}

double std_tol()
{
    return sqrt( SQR(lm_control_double.ftol) +
                 SQR(lm_control_double.gtol) +
                 SQR(lm_control_double.xtol) );
}

int check_x( int n, double *x, double *xpected, double tol )
{
    int i, ret = 0;
    for ( i=0; i<n; ++i ) {
        if ( xpected[i]==0 ) {
            if ( fabs(x[i]) > 1e-100 ) {
                printf("  x[%i] = %16.9e, xpect = 0", i, x[i] );
                ++ret;
            }
        } else {
            if ( fabs(x[i]-xpected[i]) > tol*fabs(xpected[i]) ) {
                printf("  x[%i] = %16.9e, xpect = %16.9e, relerr = %16.9e\n",
                       i, x[i], xpected[i], fabs((x[i]-xpected[i])/xpected[i]));
                ++ret;
            }
        }
    }
    return ret;
}

int check_s( int m, double *v, double spected, double tol )
{
    int j;
    double s;
    for ( j=0; j<m; ++j )
        s += v[j];
    if ( spected==0 ) {
        if ( fabs(s) > 1e-100 )
            goto failed;
    } else {
        if ( fabs(s-spected) > tol*fabs(spected) )
            goto failed;
    }
    return 0; // ok
failed:
    printf("norms deviates from expectation:\n");
    printf("           found %19.11f, expected %19.11f\n", s, spected );
    return 1;
}
