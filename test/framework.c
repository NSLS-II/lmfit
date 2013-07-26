/*
 * Library:    lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:       framework.c
 *
 * Contents:   Framework for functional tests
 *
 * Copyright:  Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 *
 * License:    see ../COPYING (FreeBSD)
 * 
 * Homepage:   apps.jcns.fz-juelich.de/lmfit
 */
 
#include "lmmin.h"
#include "framework.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#define SQR(x) ((x)*(x))

//***************************************************************************//
//  Auxiliaries
//***************************************************************************//

double std_tol()
{
    return sqrt( SQR(lm_control_double.ftol) +
                 SQR(lm_control_double.gtol) +
                 SQR(lm_control_double.xtol) );
}

//***************************************************************************//
//  Setup test cases (called by test setup functions t001, ...)
//***************************************************************************//

//! Set name of the problem (may include %-directives for parameter values).

void set_name( setup_typ* setup, const char *fmt, ... )
{
    memset( setup->name, 0, NAMELEN );
    va_list ap;
    va_start( ap, fmt );
    vsnprintf( setup->name, NAMELEN, fmt, ap );
    va_end( ap );
}

//! Set dimension of the problem; allocate parameter arrays.

void set_task( setup_typ* setup, int n_par, int m_dat, ffunc_typ f )
{
    setup->n = n_par;
    setup->m = m_dat;
    setup->f = f;
    if ( !( setup->x     = (double *) malloc( n_par * sizeof(double)) ) ||
         !( setup->xpect = (double *) malloc( n_par * sizeof(double)) ) ) {
        fprintf( stderr, "allocation of x failed\n" );
        exit( -1 );
    }
}

//! Set initial parameter values.

void set_init( setup_typ* setup, ... )
{
    va_list args;
    va_start ( args, setup );
    for ( int i = 0; i < setup->n; ++i )
        setup->x[i] = va_arg ( args, double );
    va_end ( args );
}

//! Set expected parameter values.

void set_xpec( setup_typ* setup, ... )
{
    va_list args;
    va_start ( args, setup );
    for ( int i = 0; i < setup->n; ++i )
        setup->xpect[i] = va_arg ( args, double );
    va_end ( args );
}

//***************************************************************************//
//  Register tests
//***************************************************************************//

typedef struct {
    tfunc_typ t;
    int nTP;
    double* TP;
} mini_typ;

#define MTEST 1000

mini_typ Tests[MTEST];
int nTests = 0;

void register_mini( tfunc_typ t, int nTP, ... )
{
    if( nTests+1 >= MTEST ) {
        fprintf( stderr, "Too many tests: increment MTEST\n" );
        exit(1);
    }
    Tests[nTests].t = t;
    Tests[nTests].nTP = nTP;
    if ( !( Tests[nTests].TP = (double *) malloc(nTP * sizeof(double)) ) ) {
        fprintf( stderr, "allocation of TP failed\n" );
        exit( -1 );
    }
    va_list args;
    int nargs = 0;
    va_start ( args, nTP );
    for ( int iTP = 0; iTP < nTP; ++iTP )
        Tests[nTests].TP[iTP] = va_arg ( args, double );
    va_end ( args );
    ++nTests;
}

//***************************************************************************//
//  Run registered test(s)
//***************************************************************************//

//! Private low-level implementation.

void run_test( int kTest, int verbose )
{
    printf( "%3i", kTest );
    mini_typ *T;
    T = Tests+kTest;
    setup_typ S;
    T->t( &S, T->nTP, T->TP );
    printf( " ""%-20s """, S.name );

    int n_par = S.n, m_dat = S.m;
    double tol = std_tol();

    int i, j, failed=0, badx0=0, badx1=0;
    double rel, errx0=0, errx1=0;
    static double *v;
    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    lm_princon_struct princon = lm_princon_std;
    struct timespec tim = { (time_t)0, (long)0 };
    if ( verbose ) {
        princon.form  = 1;
        princon.flags = 3;
        control.verbosity = 1;
    } else {
        princon.form  = 1;
        princon.flags = 0;
    }

    if ( !( v = (double *) realloc( v, m_dat * sizeof(double)) ) ) {
        fprintf( stderr, "allocation of v failed\n" );
        exit( -1 );
    }

    printf( "%2i %3i", n_par, m_dat );

    clock_settime( CLOCK_PROCESS_CPUTIME_ID, &tim );
    lmmin( n_par, S.x, m_dat, T->TP,
           S.f, lm_printout_std, &control, &princon, &status );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tim );
    printf( " %3i %8.4f", status.nfev, tim.tv_sec + tim.tv_nsec*1e-9 );
    if ( status.info >= 4 )
        failed = 1;

    // check fitted parameters
    for ( i=0; i<n_par; ++i ) {
        if ( S.xpect[i]==0 ) {
            if( fabs(S.x[i]) > 1e-100 ) {
                ++badx0;
                if ( fabs(S.x[i]) > errx0 )
                    errx0 = fabs(S.x[i]);
            }
        } else {
            rel = fabs((S.x[i]-S.xpect[i])/S.xpect[i]);
            if( rel > tol ) {
                ++badx1;
                if( rel > errx1 )
                    errx1 = rel;
            }
        }
    }
    if ( badx0 || badx1 )
        failed = 1;

    printf ( " %s %2i  %i %i  %8.2e %8.2e\n",
             (failed ? "failed" : "passed"), status.info,
             badx0, badx1, errx0, errx1 );
}

//! High-level API to run all tests.

void run_all()
{
     for ( int k=0; k<nTests; ++k )
         run_test( k, 0 );
}

//! High-level API to run one test.

void run_one( int kTest )
{
    if( kTest<0 || kTest>=nTests ) {
        fprintf( stderr, "invalid test number\n" );
        exit(-1);
    }
    run_test( kTest, 1 );
}
