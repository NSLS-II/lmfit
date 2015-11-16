/*
 * Library:  lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:     demo/curve1.c
 *
 * Contents: Example for curve fitting with lmcurve():
 *           fit a data set y(x) by a curve f(x;p).
 *
 * Note:     Any modification of this example should be copied to
 *           the manual page source lmcurve.pod and to the wiki.
 *
 * Author:   Joachim Wuttke <j.wuttke@fz-juelich.de> 2004-2013
 * 
 * Licence:  see ../COPYING (FreeBSD)
 * 
 * Homepage: apps.jcns.fz-juelich.de/lmfit
 */
 
#include "lmcurve.h"
#include <stdio.h>

/* model function: a parabola */

double f( double t, const double *p )
{
    return p[0] + p[1]*t + p[2]*t*t;
}

void lm_qrfac(int m, int n, double *a, int *ipvt,
              double *rdiag, double *acnorm, double *wa);

void print_matrix(int m, int n, double *a)
{
    int i,j;
    for (i=0; i<m; ++i) {
        for (j=0; j<n; ++j) {
            printf( "%8.4g ", a[j*m+i] );
        }
        printf ("\n");
    }
}

int main()
{
    double A[9] = {     0.0344,    0.4387,    0.3816,
                        0.7655,    0.7952,    0.1869,
                        0.4898,    0.4456,    0.6463 };
    int ipvt[3];
    double rdiag[3], acnorm[3], wa[3];

    printf( "Input matrix A:\n" );
    print_matrix(3, 3, A);
    
    lm_qrfac(3, 3, A, ipvt, rdiag, acnorm, wa);

    printf( "Ouput matrix A:\n" );
    print_matrix(3, 3, A);
    
    printf( "Ouput vector rdiag:\n" );
    print_matrix(1, 3, rdiag);
    
    printf( "Ouput vector acnorm:\n" );
    print_matrix(1, 3, acnorm);
}
