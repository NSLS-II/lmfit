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


int t004( int nTP, const double* TP )
{
    assert( nTP==1 );
    const int n=2, m=3;
    double xinit[2] = { 1., 1. };
    double xpect[2] = { TP[0], 2*TP[0] };
    double spect    = 0;
    return run_mini( n, m, f004, (void*)TP, xinit, xpect, spect, std_tol() );
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


int t005( int nTP, const double* TP )
{
    assert( nTP==0 );
    const int n=2, m=3;
    double xinit[2] = { 1., 1. };
    double xpect[2] = { 3., 0.5 };
    double spect    = 0;
    return run_mini( n, m, f005, (void*)TP, xinit, xpect, spect, std_tol() );
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
