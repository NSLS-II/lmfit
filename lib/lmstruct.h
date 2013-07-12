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
    double ftol;      /* Relative error desired in the sum of squares.
                         Termination occurs when both the actual and
                         predicted relative reductions in the sum of squares
                         are at most ftol. */
    double xtol;      /* Relative error between last two approximations.
                         Termination occurs when the relative error between
                         two consecutive iterates is at most xtol. */
    double gtol;      /* Orthogonality desired between fvec and its derivs.
                         Termination occurs when the cosine of the angle
                         between fvec and any column of the Jacobian is at
                         most gtol in absolute value. */
    double epsilon;   /* Step used to calculate the Jacobian, should be
                         slightly larger than the relative error in the
                         user-supplied functions. */
    double stepbound; /* Uused in determining the initial step bound. This
                         bound is set to the product of stepbound and the
                         Euclidean norm of diag*x if nonzero, or else to
                         stepbound itself. In most cases stepbound should lie
                         in the interval (0.1,100.0). Generally, the value
                         100.0 is recommended. */
    int maxcall;      /* Maximum number of minimization steps. Each step
                         requires at least n+1 function evaluations. */
    int scale_diag;   /* If 1, the variables will be rescaled internally.
                         Recommended value is 1. */
    int pivot;        /* If 1, use pivoting in QR factorization.
                         Recommended value is 1. */
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
    int info;         /* Status indicator. Nonnegative values are used as index
                         for the message text lm_infmsg, set in lmmin.c. */
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
