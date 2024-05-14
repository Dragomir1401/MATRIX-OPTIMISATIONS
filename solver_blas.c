#include "cblas.h"
#include "utils.h"

/*
 * C = ((At x B) + (B x A)) x Bt
 */
double *my_solver(int N, double *A, double *B)
{
    printf("BLAS SOLVER\n");
    double *B_copy = calloc(N * N, sizeof(double));
    double *B_copy1 = calloc(N * N, sizeof(double));
    double *C = calloc(N * N, sizeof(double));

    // Copy B to B_copy
    cblas_dcopy(N * N, B, 1, B_copy, 1);
    cblas_dcopy(N * N, B, 1, B_copy1, 1);

    // Calculate At x B (At is lower triangular)
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, 1.0, A, N, B, N);

    // Calculate B x A (A is upper triangular)
    cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1.0, A, N, B_copy, N);

    // B_copy = B x A
    // B = At x B

    // Calculate At x B + B x A
    cblas_daxpy(N * N, 1.0, B_copy, 1, B, 1);

    // B = At x B + B x A

    // Calculate ((At x B) + (B x A)) x Bt
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, B, N, B_copy1, N, 0.0, C, N);

    // Free memory
    free(B_copy);
    free(B_copy1);

    return C;
}