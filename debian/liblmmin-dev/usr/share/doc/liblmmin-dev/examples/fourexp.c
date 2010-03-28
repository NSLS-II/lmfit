/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     fourexp.c
 *
 * Contents: Fitting an t-y data set with a sum of four exponentials,
 *           as in Transtrum/Machta/Sethna PRL 104, 060201 (2010), Table I.
 *
 *           This is work in progress, a meaning ful comparison with
 *           Transtrum et al. is not yet possible.
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 
#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* model function: a sum of four exponentials */

double f( double t, const double *p )
{
    return
        exp(p[0]-exp(p[4])*t) + 
        exp(p[1]-exp(p[5])*t) + 
        exp(p[2]-exp(p[6])*t) + 
        exp(p[3]-exp(p[7])*t);
}


int main()
{
    /* initialize random number generator */
    
    srand(1);

    /* loop over random draws */

    int iter;
    for( iter=0; iter<100; ++iter ){

        /* parameter vector */

        int n_par = 8;
        double ini[8], par[8];
        int j;
        for( j=0; j<n_par; ++j ){
            ini[j] = 10. * rand()/RAND_MAX; /* to generate data */
            par[j] = 10. * rand()/RAND_MAX; /* start values for fit */
        }

        /* data pairs */

        int m_dat = 100; // number of data pairs
        double t[100], y[100];
        int i;
        for( i=0; i<m_dat; ++i ){
            /* time points chosen logarithmically spaced between 0.01 and 100 */
            t[i] = exp( log(0.01)+log(100/0.01)*i/(m_dat-1) );
            y[i] = f( t[i], ini );
        }

        /* auxiliary parameters */

        lm_status_struct status;
        lm_control_struct control = lm_control_double;
        control.maxcall = 8000;
        control.printflags = 0;

        /* perform the fit */

        lmcurve_fit( n_par, par, m_dat, t, y, f, &control, &status );

        /* print results */

        printf( "%2d %5d %18.11f\n", status.info, status.nfev, status.fnorm );

        /*
        printf( "status after %d function evaluations:\n  %s\n",
                status.nfev, lm_infmsg[status.info] );

        printf("initial and fitted parameters:\n");
        for ( i = 0; i < n_par; ++i)
            printf("  par%i: %18.11f -> %18.11f\n", i, ini[i], par[i]);
        printf("obtained norm:\n  %12g\n", status.fnorm );

        printf("fitting data as follows:\n");
        for ( i = 0; i < m_dat; ++i)
            printf( "  t[%2d]=%12g y=%12g fit=%12g residue=%12g\n",
                    i, t[i], y[i], f(t[i],par), y[i] - f(t[i],par) );
        */
    }

    return 0;
}
