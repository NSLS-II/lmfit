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
                      const lm_princon_struct *princon,
                      int iflag, int iter, int nfev );

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm( int, const double * );

/* The actual minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              const lm_princon_struct *princon,
                              int iflag, int iter, int nfev),
            const lm_control_struct *control,
            const lm_princon_struct *princon,
            lm_status_struct *status );


#ifdef __cplusplus
}
#endif

#endif /* LMMIN_H */
