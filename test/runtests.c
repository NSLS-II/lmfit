/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      runtests.c
 *
 * Contents:  Run test suite.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 *
 * License:   see ../COPYING (FreeBSD)
 * 
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void run_all();
void run_one(int);

void testsuite1();

const char* Helptext =
    "Usage:\n"
    "  runtests -h  # print help\n"
    "  runtests -a  # run all tests\n"
    "  runtests <k> # run test number k in verbose mode\n";

int main( int argc, char **argv )
{
    testsuite1(); // register tests

    if        ( argc==1 ) {
        printf( Helptext );
    } else if ( argc==2 ) {
        if      ( !strcmp(argv[1], "-h" ) ) {
            printf( Helptext ); 
            printf( "\n" );
            printf(
                "Output columns for runtests -a:\n"
                "  - serial number\n"
                "  - name of test\n"
                "  - number of parameters\n"
                "  - number of functions to be minimized\n"
                "  - duration of fit in sec\n"
                "  - return status of fit\n"
                "  - number of function calls in fit\n"
                "  - number of parameters deviating from expected value 0\n"
                "  - number of parameters deviating from expected non-0 value\n"
                "  - norm deviating from expected value 0\n"
                "  - norm deviating from expected non-0 value\n"
                "  - maximum parameter deviation from expected value 0\n"
                "  - maximum relative parameter deviation from expected non-0 value\n"
                "  - norm deviation from expected value 0\n"
                "  - norm deviation from expected non-0 value\n" );
        } else if ( !strcmp(argv[1], "-a" ) ) {
            run_all();
        } else {
            run_one( atoi(argv[1]) );
        }
    } else {
        printf( "too many arguments\n" );
        exit(-1);
    }
}
