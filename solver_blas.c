#include "utils.h"
#include "cblas.h"

// C = (At × B + B × A) × Bt

double* my_solver(int N, double *A, double *B) {
    printf("BLAS SOLVER\n");

    double *C = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));
    double *AtxB_plus_BxA = calloc(N * N, sizeof(double));

    // Calculate A^T * B where A^T is lower triangular
    // Here we use a temporary storage to hold the result of A^T * B
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, A, N, B, N, 0.0, AtxB, N);

    // Calculate B * A where A is upper triangular
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, B, N, A, N, 0.0, BxA, N);

    // Sum AtxB and BxA
    for (int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }

    // Calculate (AtxB + BxA) * B^T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, AtxB_plus_BxA, N, B, N, 0.0, C, N);

    free(AtxB);
    free(BxA);
    free(AtxB_plus_BxA);

    return C;
}
