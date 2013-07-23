int run_mini( int n_par, int m_dat,
              void (*evaluate)( const double *x, int m, const void *data,
                                double *v, int *info ),
              double *x, double* xpect, double spect, double tol );

double std_tol();
