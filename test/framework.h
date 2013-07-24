typedef int (*functyp_mini)( int, const double* );

void register_mini( const char* name, functyp_mini t, int nTP, ... );

void run_all();

int run_mini( int n_par, int m_dat,
              void (*evaluate)( const double *x, int m, const void *data,
                                double *v, int *info ),
              void* TP, double *x, double* xpect, double spect, double tol );

double std_tol();
