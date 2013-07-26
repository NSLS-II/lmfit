/*
 * Library:    lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:       testsuite1.c
 *
 * Contents:   Standard tests for minimization software.
 *
 * References: [MoGH81] JJ Moré, BS Garbow, KE Hillstrom: ACM TOMS 7, 17 (1981)
 *
 * Copyright:  Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 *
 * License:    see ../COPYING (FreeBSD)
 * 
 * Homepage:   apps.jcns.fz-juelich.de/lmfit
 */
 
#include <math.h>
#include <assert.h>
#include "framework.h"

#define SQR(x) ((x)*(x))

//  ==========================================================================
//  MoGH81 (4) Brown badly scaled function

void f004( const double *x, int m, const void *TP, double *v, int *info )
{
    assert( m== 3 );
    v[0] = x[0] - ((double*)TP)[0];
    v[1] = x[1] - 2*((double*)TP)[0];
    v[2] = x[0]*x[1] - 2;
}


void t004( setup_typ *S, int nTP, const double* TP )
{
    assert( nTP==1 );
    S->f = f004;
    set_dim( S, 2, 3 );
    set_vec( S->x,     S->n, 1., 1. );
    set_vec( S->xpect, S->n, TP[0], 2*TP[0] );
    S->spect = 0;
}

//  ==========================================================================
//  MoGH81 (5) Beale function

void f005( const double *x, int m, const void *TP, double *v, int *info )
{
    assert( m== 3 );
    v[0] = 1.5   - x[0]*(1-x[1]);
    v[1] = 2.25  - x[0]*(1-SQR(x[1]));
    v[2] = 2.625 - x[0]*(1-pow(x[1],2));
}


void t005( setup_typ *S, int nTP, const double* TP )
{
    assert( nTP==0 );
    S->f = f005;
    set_dim( S, 2, 3 );
    set_vec( S->x,     S->n, 1., 1. );
    set_vec( S->xpect, S->n, 3., 0.5 );
    S->spect = 0;
}

//  ==========================================================================
//  register examples

int testsuite1()
{
    register_mini( "MoGH81#04", t004, 1, 1.e3 );
    register_mini( "MoGH81#04", t004, 1, 1.e6 );
    register_mini( "MoGH81#04", t004, 1, 1.e12 );
    register_mini( "MoGH81#04", t004, 1, 1.e18 );

    register_mini( "MoGH81#05", t005, 0 );
}
