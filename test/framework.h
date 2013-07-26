typedef void (*ffunc_typ)( const double*x, int, const void*, double*, int* );

#define NAMELEN 20

typedef struct {
    char name[NAMELEN];
    int m;
    int n;
    ffunc_typ f;
    double *x;
    double *xpect;
    double spect;
} setup_typ;

void set_name( setup_typ* setup, const char *fmt, ... );
void set_task( setup_typ* setup, int n_par, int m_dat, ffunc_typ f );
void set_init( setup_typ* setup, ... );
void set_xpec( setup_typ* setup, ... );

typedef void (*tfunc_typ)( setup_typ*, int, const double* );

void register_mini( tfunc_typ t, int nTP, ... );

void run_all();
void run_one( int kTest );

int run_mini( int n_par, int m_dat,
              void (*evaluate)( const double *x, int m, const void *data,
                                double *v, int *info ),
              void* TP, double *x, double* xpect, double spect, double tol );

double std_tol();
