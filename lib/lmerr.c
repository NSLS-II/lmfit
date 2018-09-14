/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmerr.c
 *
 * Contents:  Error estimation.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2018-)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */

#include "lmerr.h"
#include "lmmin.h" /* for fnorm */
#include <math.h>
#include <stdlib.h>
#include <float.h>

#define SQR(x)   (x)*(x)
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))
/* machine-dependent constants from float.h */
#define LM_MACHEP     DBL_EPSILON   /* resolution of arithmetic */
#define LM_DWARF      DBL_MIN       /* smallest nonzero number */

/*****************************************************************************/
/*  lmerr (diagonal error estimates)                                         */
/*****************************************************************************/

void lmerr(
    const int n, double *const x, double *const dx,
    const int m, const double* y, const void *const data,
    void (*const evaluate)(
        const double *const par, const int m_dat, const void *const data,
        double *const fvec, int *const userbreak),
    const lm_control_struct *const C, int* failure)
{
    int j, i;
    double temp, fnorm, step, df, sum;
    double eps = sqrt(MAX(C->epsilon, LM_MACHEP)); /* for forward differences */
    *failure = 0;

/***  Check input parameters for errors.  ***/

    if ( n < 0 ) {
        *failure = 10; /* invalid parameter */
        return;
    }
    if (m < n) {
        *failure = 10;
        return;
    }

/***  Allocate work space.  ***/

    /* Allocate total workspace with just one system call */
    char *ws;
    if ( ( ws = malloc((2*m)*sizeof(double)) ) == NULL ) {
        *failure = 9;
        return;
    }

    /* Assign workspace segments. */
    char *pws = ws;
    double *fvec = (double*) pws; pws += m*sizeof(double)/sizeof(char);
    double *wf   = (double*) pws; pws += m*sizeof(double)/sizeof(char);

/***  Calculate the Jacobian.  ***/

    (*evaluate)(x, m, data, fvec, failure);
    if ( *failure )
        return;
    fnorm = lm_fnorm(m, fvec, y);

    if( fnorm <= LM_DWARF ) {
        for (j = 0; j < n; j++)
            dx[j] = 0.;

    } else {
        for (j = 0; j < n; j++) {
            sum = 0;
            temp = x[j];
            step = MAX(eps*eps, eps * fabs(temp));
            x[j] += step; /* replace temporarily */
            (*evaluate)(x, m, data, wf, failure);
            if ( *failure )
                return;
            for (i = 0; i < m; i++) {
                df = (wf[i] - fvec[i]) / step;
                sum += SQR(df/fnorm);
            }
            dx[j] = 1 / sqrt(sum);
            x[j] = temp; /* restore */
        }
    }

/***  Cleanup.  ***/
    free(ws);
}
