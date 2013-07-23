int run_mini( int n_par, int m_dat,
              void (*evaluate)( const double *x, int m, const void *data,
                                double *v, int *info ),
              double *x, double *v );

double std_tol();

int check_x( int n, double *x, double *xpected, double tol );

int check_s( int m, double *v, double spected, double tol );
