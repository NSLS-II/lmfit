/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmcurve.h
 *
 * Contents: Simplified interface for one-dimensional curve fitting
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  CC-BY-SA (Creative Commons Attribution Share-Alike)
 */
 
#include<lmmin.h>

#ifndef LMCURVE_H
#define LMCURVE_H

#ifdef __cplusplus
extern "C" {
#endif

void lmcurve_fit( int m_dat, int n_par, double *par,
                  double *t, double *y, double (*f)( double t, double *par ),
                  lm_control_type *control );

#ifdef __cplusplus
}
#endif

#endif /* LMCURVE_H */
