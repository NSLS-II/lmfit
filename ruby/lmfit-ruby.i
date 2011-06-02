%module lmfit

%{
#include <lmmin.h>
#include <lmcurve.h>
%}

%include "cpointer.i"
%include "carrays.i"

%pointer_functions(unsigned short, usp)
%array_functions(unsigned char, uca)

void lmcurve_fit( int, double*, int, const double*, const double*,
                  double (*f)( double, const double *),
                  lm_control_struct*, lm_status_struct* );
