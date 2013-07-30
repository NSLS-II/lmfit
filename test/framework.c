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
    // Call Tests[kTest].t() to obtain test setup S for given parameters TP.
    mini_typ *T;
    T = Tests+kTest;
    setup_typ S;
    T->t( &S, T->nTP, T->TP );
    printf( "%3i %-20s %2i %3i", kTest, S.name, S.n, S.m );

    // Prepare for the minimization.
    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    lm_princon_struct princon = lm_princon_std;
    if ( verbose ) {
        princon.form  = 1;
        princon.flags = 3;
        control.verbosity = 1;
    } else {
        princon.form  = 1;
        princon.flags = 0;
    }
    struct timespec tim = { (time_t)0, (long)0 };

    // Run the minization.
    if (verbose )
        printf( ":\n" );
    clock_settime( CLOCK_PROCESS_CPUTIME_ID, &tim );
    lmmin( S.n, S.x, S.m, T->TP,
           S.f, lm_printout_std, &control, &princon, &status );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tim );
    if (verbose )
        printf( "------------------------------------> " );

    // Check fitted parameters.
    int badx=0;
    double errx=0, rel;
    for ( int i=0; i<S.n; ++i ) {
        if ( fabs(S.xpect[i]) < 1e-90 ?
             ( rel = fabs(S.x[i]) ) > 1e-99 :
             ( rel = fabs((S.x[i]-S.xpect[i])/S.xpect[i]) ) > 1e-14 ) {
            ++badx;
            if( rel > errx )
                errx = rel;
        }
    }
    const char *result;
    if      ( status.info >= 4 )
        result = "failed";
    else if ( badx )
        result = "doubt";
    else
        result = "passed";
        
    // Print outcome.
    printf( " %8.4f %1i %3i %6s %2i %8.2e\n", tim.tv_sec + tim.tv_nsec*1e-9,
            status.info, status.nfev, result, badx, errx );
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
