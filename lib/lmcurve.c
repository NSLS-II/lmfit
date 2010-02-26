/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmcurve.c
 *
 * Contents: Simplified wrapper for one-dimensional curve fitting,
 *           using Levenberg-Marquardt least-squares minimization.
 *
 * Usage:    see application sample demo/curve1.c
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Creative Commons Attribution Share Alike.
 */
 

#include "lmmin.h"


typedef struct {
    double *t;
    double *y;
    double (*f) (double t, double *par);
} lmcurve_data_type;


void lmcurve_evaluate( double *par, int m_dat, double *fvec,
		       void *data, int *info )
{
    int i;
    for ( i = 0; i < m_dat; i++ )
	fvec[i] =
            ((lmcurve_data_type*)data)->y[i] -
            ((lmcurve_data_type*)data)->f(
                ((lmcurve_data_type*)data)->t[i], par );
    *info = *info;		/* to prevent a 'unused variable' warning */
}


void lmcurve_fit( int m_dat, int n_par, double *par,
                  double *t, double *y, double (*f)( double t, double *par ),
                  lm_control_type *control )
{
    lmcurve_data_type data;
    data.t = t;
    data.y = y;
    data.f = f;
    lm_minimize( m_dat, n_par, par, lmcurve_evaluate, 0,
                 (void*) &data, control );
}

