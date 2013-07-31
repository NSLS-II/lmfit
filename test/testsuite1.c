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

void f004( const double *x, int m, const void *TP, double *v, int *usrbrk )
{
    assert( m==3 );
    v[0] = x[0] - ((double*)TP)[0];
    v[1] = x[1] - 2/((double*)TP)[0];
    v[2] = x[0]*x[1] - 2;
}

void t004( setup_typ *S, int nTP, const double* TP )
{
    assert( nTP==1 );
    set_name( S, "MoGH81#04[%6.1e]", TP[0] );
    set_task( S, 2, 3, f004 );
    set_init( S, 1., 1. );
    set_xpec( S, TP[0], 2/TP[0] );
}

//  ==========================================================================
//  MoGH81 (5) Beale function

void f005( const double *x, int m, const void *TP, double *v, int *usrbrk )
{
    assert( m==3 );
    v[0] = 1.5   - x[0]*(1-x[1]);
    v[1] = 2.25  - x[0]*(1-SQR(x[1]));
    v[2] = 2.625 - x[0]*(1-pow(x[1],3));
}

void t005( setup_typ *S, int nTP, const double* TP )
{
    assert( nTP==0 );
    set_name( S, "MoGH81#05" );
    set_task( S, 2, 3, f005 );
    set_init( S, 1., 1. );
    set_xpec( S, 3., 0.5 );
}


//  ==========================================================================
//  MoGH81 (32) linear function - full rank

void f032( const double *x, int m, const void *TP, double *v, int *usrbrk )
{
    int n = lrint( ((double*)TP)[0] );
    double s = 0;
    for ( int j=0; j<n; ++j )
        s += x[j];
    s = -2*s/m-1;
    for ( int i=0; i<n; ++i )
        v[i] = x[i] + s;
    for ( int i=n; i<m; ++i )
        v[i] = s;
}

void t032( setup_typ *S, int nTP, const double* TP )
{
    assert( nTP==2 );
    int n = lrint( TP[0] );
    int m = lrint( TP[1] );
    set_name( S, "MoGH81#32[%i,%i]", n, m );
    set_task( S, n, m, f032 );
    for ( int j=0; j<n; ++j ) {
        S->x[j] = 1;
        S->xpect[j] = -1;
    }
}


//  ==========================================================================
//  register examples

int testsuite1()
{
    register_mini( t004, 1, 1.e3 );
    register_mini( t004, 1, 1.e6 );
    register_mini( t004, 1, 1.e7 );
    register_mini( t004, 1, 1.e8 );

    register_mini( t005, 0 );

    register_mini( t032, 2, 4., 7. );
}
