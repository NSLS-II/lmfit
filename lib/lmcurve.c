/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve.c
 *
 * Contents:  Implements lmcurve, a simplified API for curve fitting
 *            using the generic Levenberg-Marquardt routine lmmin.
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

/******************************************************************************/
/*  lmcurve                                                                   */
/******************************************************************************/

typedef struct {
    const double* t;
    const double* y;
    double (*f)(const double t, const double* par);
} lmcurve_data_struct;

void lmcurve_evaluate(
    const double* par, const int m_dat, const void* data, double* fvec,
    int* info)
{
    lmcurve_data_struct* D = (lmcurve_data_struct*)data;
    int i;
    for (i = 0; i < m_dat; i++)
        fvec[i] = D->y[i] - D->f(D->t[i], par);
}

void lmcurve(
    const int n_par, double* par, const int m_dat,
    const double* t, const double* y,
    double (*f)(const double t, const double* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_data_struct data = { t, y, f };

    lmmin(n_par, par, m_dat, (const void*)&data, lmcurve_evaluate,
          control, status);
}

/******************************************************************************/
/*  lmcurve_tyd                                                               */
/******************************************************************************/

typedef struct {
    const double* t;
    const double* y;
    const double* dy;
    double (*f)(const double t, const double* par);
} lmcurve_tyd_data_struct;

void lmcurve_tyd_evaluate(
    const double* par, const int m_dat, const void* data, double* fvec,
    int* info)
{
    lmcurve_tyd_data_struct* D = (lmcurve_tyd_data_struct*)data;
    int i;
    for (i = 0; i < m_dat; i++)
        fvec[i] = ( D->y[i] - D->f(D->t[i], par) ) / D->dy[i];
}

void lmcurve_tyd(
    const int n_par, double* par, const int m_dat,
    const double* t, const double* y, const double* dy,
    double (*f)(const double t, const double* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_tyd_data_struct data = { t, y, dy, f };

    lmmin(n_par, par, m_dat, (const void*)&data, lmcurve_tyd_evaluate,
          control, status);
}
