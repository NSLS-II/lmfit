/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     curve1.c
 *
 * Contents: Example for one-dimensional curve fitting,
 *           using the simplified interface from lmcurve.h.
 *
 * Author:   Joachim Wuttke 2004-10
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Creative Commons Attribution Share Alike.
 */
 
#include "lmcurve.h"
#include <stdio.h>

double f(double t, double *p)
{
    return (p[0] * t + (1 - p[0] + p[1] + p[2]) * t * t) /
	(1 + p[1] * t + p[2] * t * t);
}

int main()
{
    // data and pameter arrays:

    int m_dat = 15;
    int n_par = 3;
    int i;

    double t[15] = { .07, .13, .19, .26, .32, .38, .44, .51,
	.57, .63, .69, .76, .82, .88, .94
    };
    double y[15] = { .24, .35, .43, .49, .55, .61, .66, .71,
	.75, .79, .83, .87, .90, .94, .97
    };
    double par[3] = { 1., 1., 1. }; // use any starting value, but not { 0,0,0 }

    // auxiliary parameter records:

    lm_control_struct control = lm_control_double;
    lm_status_struct status;

    // maximum noise:

    control.printflags = 7;

    // perform the fit:

    lmcurve_fit( m_dat, n_par, par, t, y, f, &control, &status );

    // print results:

    printf( "status after %d evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.info] );

    printf("obtained parameters:\n");
    for (i = 0; i < n_par; ++i)
	printf("  par[%i] = %12g\n", i, par[i]);
    printf("obtained norm:\n  %12g\n", status.fnorm );

    printf("fitting data as follows:\n");
    for (i = 0; i < m_dat; ++i)
        printf( "  t[%2d]=%12g y=%12g fit=%12g residue=%12g\n",
                i, t[i], y[i], f(t[i],par), y[i] - f(t[i],par) );
    return 0;
}
