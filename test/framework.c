/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     framework.c
 *
 * Contents: Framework for functional tests.
 *
 * Author:   Joachim Wuttke 2013
 * 
 * Homepage: apps.jcns.fz-juelich.de/lmfit
 *
 * Licence:  see ../COPYING (FreeBSD)
 */
 
#include "lmmin.h"
#include "framework.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define SQR(x) ((x)*(x))

typedef struct {
    const char* name;
    functyp_mini t;
    int nTP;
    double* TP;
} mini_typ;

#define MTEST 1000

mini_typ Tests[MTEST];
int nTests = 0;

void register_mini( const char* name, functyp_mini t, int nTP, ... )
{
    if( nTests+1 >= MTEST ) {
        fprintf( stderr, "Too many tests: increment MTEST\n" );
        exit(1);
    }
    Tests[nTests].name = name;
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

void run_test( int m )
{
    printf( "TEST %3i", m );
    mini_typ *T;
    T = Tests+m;
    printf( " ""%-12s""", T->name );
    T->t( T->nTP, T->TP );
}

void run_all()
{
     for ( int m=0; m<nTests; ++m )
         run_test( m );
}


int run_mini( int n_par, int m_dat,
              void (*evaluate)( const double *x, int m, const void *data,
                                double *v, int *info ),
              void* TP, double *x, double* xpect, double spect, double tol )
{
    int i, j, failed = 0, badx = 0, bads = 0;
    double s;
    double *v;
    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    lm_princon_struct princon = lm_princon_std;
    struct timespec tim = { (time_t)0, (long)0 };
    princon.form  = 1;
    princon.flags = 0;

    if ( !( v = (double *) malloc(m_dat * sizeof(double)) ) ) {
        fprintf( stderr, "allocation of v failed\n" );
        exit( -1 );
    }

    clock_settime( CLOCK_PROCESS_CPUTIME_ID, &tim );
    lmmin( n_par, x, m_dat, TP,
           evaluate, NULL /*lm_printout_std*/, &control, &princon, &status );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tim );
    printf( " %3i %8.4f", status.nfev, tim.tv_sec + tim.tv_nsec*1e-9 );
    if ( status.info >= 4 )
        failed = 1;

    // check fitted parameters
    for ( i=0; i<n_par; ++i ) {
        if ( xpect[i]==0 ?
             fabs(x[i]) > 1e-100 :
             fabs(x[i]-xpect[i]) > tol*fabs(xpect[i]) )
            ++badx;
    }
    if ( badx )
        failed = 1;

    // check obtained minimum
    for ( j=0; j<m_dat; ++j )
        s += v[j];
    if ( spect==0 ?
         fabs(s) > 1e-100 :
         fabs(s-spect) > tol*fabs(spect) ) {
        bads = 1;
        failed = 1;
    }

    printf ( " %s %2i %i %i\n",
             failed ? "failed" : "passed", status.info, badx, bads );

    return failed;
}

double std_tol()
{
    return sqrt( SQR(lm_control_double.ftol) +
                 SQR(lm_control_double.gtol) +
                 SQR(lm_control_double.xtol) );
}
