typedef void (*ffunc_typ)( const double*x, int, const void*, double*, int* );

#define NAMELEN 20

typedef struct {
    char name[NAMELEN];
    int m;
    int n;
    ffunc_typ f;
    double *x;
    double *xpect;
} setup_typ;

void set_name( setup_typ* setup, const char *fmt, ... );
void set_task( setup_typ* setup, int n_par, int m_dat, ffunc_typ f );
void set_init( setup_typ* setup, ... );
void set_xpec( setup_typ* setup, ... );

typedef void (*tfunc_typ)( setup_typ*, int, const double* );

void register_mini( tfunc_typ t, int nTP, ... );

double std_tol();
