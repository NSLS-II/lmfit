/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     powell.c
 *
 * Contents: A problem from Powell 1970, according to Madsen et al. 2004,
 *           Example 3.2
 *
 *           Minimize the norm of a 2-dimensional vectorial function
 *           of a 2-dimensional parameter vector:
 *           F(p) = ( p0, 10*p0/(p0+0.1)+2*p1^2 ).
 *
 *           The only analytical solution is par = ( 0, 0 ).
 *           The problem is nontrivial because the Jacobian is
 *           singular at the solution. And with starting value ( 3, 1 ),
 *           exact line search converges to a point that is no solution.
 *
 *           With the present LM implementation, on a x86 slow but stable
 *           convergence is found.
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 
#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>


void evaluate_powell( const double *par, int m_dat, const void *data,
                      double *fvec, int *info )
{
    fvec[0] = par[0];
    fvec[1] = 10*par[0]/(par[0]+0.1) + 2*par[1]*par[1];
}


int main()
{
    /* parameter vector */

    int n_par = 2; // number of parameters in model function f
    double par[2] = { 3, 1 }; // arbitrary starting value

    /* data points */

    int m_dat = 2;

    /* auxiliary parameters */

    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    control.printflags = 3; // monitor status (+1) and parameters (+2)

    /* perform the fit */

    printf( "Fitting:\n" );
    lmmin( n_par, par, m_dat, 0,
           evaluate_powell, &control, &status, lm_printout_std );

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
