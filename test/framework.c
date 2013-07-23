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
              double *x, double* xpect, double spect, double tol )
{
    int i, ret = 0;
    int j;
    double s;
    double *v;
    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    lm_princon_struct princon = lm_princon_std;
    princon.form  = 1;
    princon.flags = 0;

    if ( ( v = (double *) malloc(m_dat * sizeof(double))) == NULL ) {
        fprintf( stderr, "allocation of v failed\n" );
        exit( -1 );
    }
    lmmin( n_par, x, m_dat, NULL,
           evaluate, NULL /*lm_printout_std*/, &control, &princon, &status );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.info] );
    if ( status.info >= 4 )
        return status.info;

    // check fitted parameters
    for ( i=0; i<n_par; ++i ) {
        if ( xpect[i]==0 ) {
            if ( fabs(x[i]) > 1e-100 ) {
                printf("  x[%i] = %16.9e, xpect = 0", i, x[i] );
                ++ret;
            }
        } else {
            if ( fabs(x[i]-xpect[i]) > tol*fabs(xpect[i]) ) {
                printf("  x[%i] = %16.9e, xpect = %16.9e, relerr = %16.9e\n",
                       i, x[i], xpect[i], fabs((x[i]-xpect[i])/xpect[i]));
                ++ret;
            }
        }
    }
    if ( ret )
        return ret;

    // check obtained minimum
    for ( j=0; j<m_dat; ++j )
        s += v[j];
    if ( spect==0 ) {
        if ( fabs(s) > 1e-100 )
            goto failed;
    } else {
        if ( fabs(s-spect) > tol*fabs(spect) )
            goto failed;
    }
    return 0; // ok
failed:
    printf("norms deviates from expectation:\n");
    printf("           found %19.11f, expect %19.11f\n", s, spect );
    return 1;
}

double std_tol()
{
    return sqrt( SQR(lm_control_double.ftol) +
                 SQR(lm_control_double.gtol) +
                 SQR(lm_control_double.xtol) );
}
