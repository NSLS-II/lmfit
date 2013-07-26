/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmmin.h
 *
 * Contents:  Declarations for Levenberg-Marquardt minimization.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 * 
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */

#ifndef LMMIN_H
#define LMMIN_H

#include "lmstruct.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*lm_eval_ftyp) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info);
typedef void (*lm_prin_ftyp) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              const lm_princon_struct *princon,
                              int iflag, int iter, int nfev);

/* Levenberg-Marquardt minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            lm_eval_ftyp evaluate,
            lm_prin_ftyp printout,
            const lm_control_struct *control,
            const lm_princon_struct *princon,
            lm_status_struct *status );
/*
 *   This routine contains the core algorithm of our library.
 *
 *   It minimizes the sum of the squares of m nonlinear functions
 *   in n variables by a modified Levenberg-Marquardt algorithm.
 *   The function evaluation is done by the user-provided routine 'evaluate'.
 *   The Jacobian is then calculated by a forward-difference approximation.
 *
 *   Parameters:
 *
 *      n is the number of variables (INPUT, positive integer).
 *
 *      x is the solution vector (INPUT/OUTPUT, array of length n).
 *        On input it must be set to an estimated solution.
 *        On output it yields the final estimate of the solution.
 *
 *      m is the number of functions to be minimized (INPUT, positive integer).
 *        It must fulfill m>=n.
 *
 *      data is a pointer that is ignored by lmmin; it is however forwarded
 *        to the user-supplied functions evaluate and printout.
 *        In a typical application, it contains experimental data to be fitted.
 *
 *      evaluate is a user-supplied function that calculates the m functions.
 *        Parameters:
 *          n, x, m, data as above.
 *          fvec is an array of length m; on OUTPUT, it must contain the
 *            m function values for the parameter vector x.
 *          info is an integer pointer. When *info is set to a negative value,
 *            lmmin will terminate.
 *
 *      printout is a user-supplied function that informs about fit progress.
 *        Call with printout=NULL if no printout is desired.
 *        Call with printout=lm_printout_std to use the default implementation.
 *
 *      control contains INPUT variables that control the fit algorithm,
 *        as declared and explained in lmstruct.h
 *
 *      princon contains INPUT variables that control the monitoring,
 *        as declared and explained in lmstruct.h
 *
 *      status contains OUTPUT variables that inform about the fit result,
 *        as declared and explained in lmstruct.h
 */

/* Standard monitoring routine. */
void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      const lm_princon_struct *princon,
                      int iflag, int iter, int nfev );

/* Refined calculation of Eucledian norm. */
double lm_enorm( int, const double * );

#ifdef __cplusplus
}
#endif

#endif /* LMMIN_H */
