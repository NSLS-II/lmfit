#include "lminvert.h"
#include <stdio.h>
#include <stdlib.h>

void show_matrix(double *A, int n) {
    printf("\n");
    for (int j=0; j < n; j++) {
        for (int i=0; i < n; i++)
            printf("%12.5f ", A[j*n+i]);
        printf("\n");
    }
}

void test(int n, double* A)
{
    double* ws = (double*)calloc(n*n, sizeof(double));
    double* inv = (double*)calloc(n*n, sizeof(double));
    int*    P = (int*)calloc(n, sizeof(int));
    if (!ws || !inv)
        exit(1);
    show_matrix(A, n);
    int failure = 0;
    lm_invert(A, n, P, ws, inv, &failure);
    if (failure==0) {
        ; // ok
    } else if (failure==21) {
        printf("singular\n");
        return;
    } else if (failure==22) {
        printf("inaccurate\n");
        return;
    } else {
        printf("unexpected error %i\n", failure);
        return;
    }
    show_matrix(inv, n);
}

int main() {
    double m0[] = {1, 1.,
                   0, 1};
    test(2, m0);

    double ms[] = {1, 3.,
                   1, 3};
    test(2, ms);

    double m1[] = {25, 15, -5,
                   15, 18,  0,
                   -5,  0, 11};
    test(3, m1);

    double m2[] = {18, 22,  54,  42,
                   22, 70,  86,  62,
                   54, 86, 174, 134,
                   42, 62, 134, 106};
    test(4, m2);

    return 0;
}
