/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve.c
 *
 * Contents:  Levenberg-Marquardt curve-fitting
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 *
 * Note to programmers: Don't patch and fork, but copy and variate!
 *   If you need to compute residues differently, then please do not patch
 * lmcurve.h and lmcurve.c, but copy them, and create differently named
 * versions of lmcurve_data_struct, lmcurve_evaluate, and lmcurve of your own.
 */

#include "lmmin.h"


typedef struct {
    const double *t;
    const double *y;
    double (*f) (double t, const double *par);
} lmcurve_data_struct;


void lmcurve_evaluate( const double *par, int m_dat, const void *data,
                       double *fvec, int *info )
{
    int i;
    for ( i = 0; i < m_dat; i++ )
        fvec[i] =
            ((lmcurve_data_struct*)data)->y[i] -
            ((lmcurve_data_struct*)data)->f(
                ((lmcurve_data_struct*)data)->t[i], par );
}


void lmcurve( int n_par, double *par, int m_dat,
              const double *t, const double *y,
              double (*f)( double t, const double *par ),
              const lm_control_struct *control,
              lm_status_struct *status )
{
    lmcurve_data_struct data;
    data.t = t;
    data.y = y;
    data.f = f;

    lmmin( n_par, par, m_dat, (const void*) &data,
           lmcurve_evaluate, control, status );
}
