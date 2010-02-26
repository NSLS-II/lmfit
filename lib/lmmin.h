/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmmin.h
 *
 * Contents: Public interface to the Levenberg-Marquardt core implementation.
 *
 * Author:   Joachim Wuttke 2004-2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Creative Commons Attribution Share Alike.
 */
 
#ifndef LMMIN_H
#define LMMIN_H

#ifdef __cplusplus
extern "C" {
#endif


/** Compact high-level interface. **/

/* Collection of control and monitoring parameters. */
typedef struct {
    /* control (input) parameters */
    double ftol;      /* relative error desired in the sum of squares. */
    double xtol;      /* relative error between last two approximations. */
    double gtol;      /* orthogonality desired between fvec and its derivs. */
    double epsilon;   /* step used to calculate the jacobian. */
    double stepbound; /* initial bound to steps in the outer loop. */
    int maxcall;      /* maximum number of iterations. */
    /* monitoring (output) parameters */
    double fnorm;     /* norm of the residue vector fvec. */
    int nfev;	      /* actual number of iterations. */
    int info;	      /* status of minimization. */
} lm_control_type;

/* Recommended parameter settings. */
extern lm_control_type lm_control_double;
extern lm_control_type lm_control_float;

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm(int, double *);

/* The actual minimization. */
void lm_minimize(int m_dat, int n_par, double *par,
		 void (*evaluate) (double *par, int m_dat, double *fvec,
                                   void *data, int *info),
                 void (*printout) (int n_par, double *par, int m_dat,
                                   double *fvec, void *data, int iflag,
                                   int iter, int nfev),
		 void *data, lm_control_type * control);


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lm_minimize. */
void lm_lmdif(int m, int n, double *x, double *fvec, double ftol,
	      double xtol, double gtol, int maxfev, double epsfcn,
	      double *diag, int mode, double factor, int *info, int *nfev,
	      double *fjac, int *ipvt, double *qtf, double *wa1,
	      double *wa2, double *wa3, double *wa4,
              void (*evaluate) (double *par, int m_dat, double *fvec,
                                void *data, int *info),
              void (*printout) (int n_par, double *par, int m_dat,
                                double *fvec, void *data, int iflag,
                                int iter, int nfev),
	      void *data);

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];


/*** the following is OBSOLETE since version 3.0 ***/

/* Initialize control parameters with default values. */
void lm_initialize_control(lm_control_type * control);

/* Call-back routines for one-dimensional curve fitting. */

void lm_evaluate_default(double *par, int m_dat, double *fvec, void *data,
			 int *info);

void lm_print_default(int n_par, double *par, int m_dat, double *fvec,
		      void *data, int iflag, int iter, int nfev);

/* Record type for passing data and model function to lm_evaluate. */

typedef struct {
    double *tvec;
    double *yvec;
    double (*f) (double t, double *par);
} lm_data_type_default;


#ifdef __cplusplus
}
#endif

#endif /* LMMIN_H */
