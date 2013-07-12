/*
 * Library:  lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:     lmcurve.h
 *
 * Contents: Simplified interface for one-dimensional curve fitting
 *
 * Author:   Joachim Wuttke <j.wuttke@fz-juelich.de> 2004-2013
 * 
 * Licence:  see ../COPYING (FreeBSD)
 * 
 * Homepage: apps.jcns.fz-juelich.de/lmfit
 */
 
#ifndef LMCURVE_H
#define LMCURVE_H

#include<lmstruct.h>

#ifdef __cplusplus
extern "C" {
#endif

void lmcurve( int n_par, double *par, int m_dat,
              const double *t, const double *y,
              double (*f)( double t, const double *par ),
              const lm_control_struct *control,
              const lm_princon_struct *princon,
              lm_status_struct *status );

#ifdef __cplusplus
}
#endif

#endif /* LMCURVE_H */
