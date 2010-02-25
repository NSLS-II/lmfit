/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmfit.h
 *
 * Contents: Public interface to the Levenberg-Marquardt core implementation.
 *
 * Author:   Joachim Wuttke 2004-2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Public domain.
 */
 
#include<lmmin.h>

#ifndef LMFIT_H
#define LMFIT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Call-back routines for one-dimensional curve fitting: */

void lmfit_evaluate( double *par, int m_dat, double *fvec,
                     void *data, int *info );

void lmfit_printout( int n_par, double *par, int m_dat, double *fvec,
		     void *data, int iflag, int iter, int nfev );


/* Record type for passing data and model function to lmfit_evaluate: */

typedef struct {
    double *tvec;
    double *yvec;
    double (*f)( double t, double *par );
} lmfit_data_type;

#ifdef __cplusplus
}
#endif

#endif /* LMFIT_H */
