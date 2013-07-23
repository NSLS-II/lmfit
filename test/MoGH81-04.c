/* MoGH81 (4) Brown badly scaled function */

#include <assert.h>
#include "framework.h"

void evaluate ( const double *x, int m, const void *data,
                double *v, int *info )
{
    assert( m== 3 );
    v[0] = x[0] - 1e6;
    v[1] = x[1] - 2e-6;
    v[2] = x[0]*x[1] - 2;
}

int main()
{
    const int n=2, m=3;
    double x[2] = { 1., 1. };
    double X[2] = { 1e6, 2e-6 };
    double v[3];
    int ret;
    if( ret = run_mini( n, m, evaluate, x, v ) )
        return ret;
    if( check_x( n, x, X, std_tol() ) )
        return 20;
    if( check_s( m, v, 0., std_tol() ) )
        return 30;
    return 0;
}
