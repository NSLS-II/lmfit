/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmmin.h
 *
 * Contents: Public interface to the Levenberg-Marquardt core implementation.
 *
 * Author:   Joachim Wuttke 2004-2013
 *
 * Licence:  see ../COPYING (FreeBSD)
 * 
 * Homepage: joachimwuttke.de/lmfit
 */
 
#ifndef LMMIN_H
#define LMMIN_H

#include<lmstruct.h>

#ifdef __cplusplus
extern "C" {
#endif


/* Standard monitoring routine. */
void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev );

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm( int, const double * );

/* The actual minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) );


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lmmin. */
void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
               double xtol, double gtol, int maxfev, double epsfcn,
               double *diag, int mode, double factor, int *info, int *nfev,
               double *fjac, int *ipvt, double *qtf, double *wa1,
               double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
               int printflags, const void *data );

#ifdef __cplusplus
}
#endif

#endif /* LMMIN_H */
