/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmerr.h
 *
 * Contents:  Declarations for error estimation.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2018-)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */

#ifndef LMERR_H
#define LMERR_H
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

#include "lmstruct.h"

__BEGIN_DECLS

/* Estimation of parameter uncertainties. */
void lmerr(
    const int n_par, double* par, double* parerr,
    const int m_dat, const double* y, const void* data,
    void (*evaluate)(
        const double* par, const int m_dat, const void* data,
        double* fvec, int* userbreak),
    const lm_control_struct* control, int* failure);
/*
 *   To be called after successful least-squares minimization.
 *   Computes parerr.
 *   Note that the very concept of parameter uncertainty (= error bar)
 *   is too naive since it ignores covariances.
 *   Here, we just use the diagonal elements of the covariance matrix.
 *
 *   Parameters:
 *
 *      n_par is the number of variables (INPUT, positive integer).
 *
 *      par is the parameter vector (INPUT, array of length n) obtained by lmmin.
 *
 *      parerr will contain the error estimates (OUTPUT, array of length n).
 *
 *      m_dat is the number of functions to be minimized (INPUT, integer).
 *        It must fulfill m>=n.
 *
 *      y contains data to be fitted. Use a null pointer if there are no data.
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
 *          userbreak is an integer pointer. When *userbreak is set to a
 *            nonzero value, lmmin will terminate.
 *
 *      control contains INPUT variables that control the fit algorithm,
 *        as declared and explained in lmstruct.h
 *
 *      *failure (OUTPUT, integer) is set to 0 unless the computation failed.
 */

__END_DECLS
#endif /* LMERR_H */
