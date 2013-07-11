/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmstruct.h
 *
 * Contents: Data structures containing control and status variables.
 *
 * Author:   Joachim Wuttke 2004-2013
 *
 * Licence:  see ../COPYING (FreeBSD)
 * 
 * Homepage: joachimwuttke.de/lmfit
 */
 
#ifndef LMSTRUCT_H
#define LMSTRUCT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Collection of input parameters for fit control. */
typedef struct {
    double ftol;      /* relative error desired in the sum of squares. */
    double xtol;      /* relative error between last two approximations. */
    double gtol;      /* orthogonality desired between fvec and its derivs. */
    double epsilon;   /* step used to calculate the jacobian. */
    double stepbound; /* initial bound to steps in the outer loop. */
    int maxcall;      /* maximum number of iterations. */
    int scale_diag;   /* automatical diag rescaling? [UNDOCUMENTED] */
    int pivot;        /* use pivoting in QR factorization? [UNDOCUMENTED] */
} lm_control_struct;

/* Collection of input parameters for print control. */
typedef struct {
    int form;         /* to select one out of several forms. */
    int flags;        /* OR'ed switches that decide which info gets printed. */
    int n_maxpri;     /* -1, or max number of parameters to print. */
    int m_maxpri;     /* -1, or max number of residuals to print. */
} lm_princon_struct;

/* Collection of output parameters for status info. */
typedef struct {
    double fnorm;     /* norm of the residue vector fvec. */
    int nfev;         /* actual number of iterations. */
    int info;         /* status (index for lm_infmsg and lm_shortmsg). */
} lm_status_struct;

/* Preset (and recommended) control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;

extern const lm_princon_struct lm_princon_std;

/* Preset message texts. */

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];


#ifdef __cplusplus
}
#endif

#endif /* LMSTRUCT_H */
