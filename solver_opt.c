/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"

// C = (At × B + B × A) × Bt

/*
 * Add your optimized implementation here
 */
void transpose(double *src, double *dst, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            dst[j * N + i] = src[i * N + j];
        }
    }
}

double *my_solver(int N, double *A, double *B) {
    printf("OPT SOLVER\n");

    double *At = calloc(N * N, sizeof(double));
    double *Bt = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));
    double *AtxB_plus_BxA = calloc(N * N, sizeof(double));
    double *C = calloc(N * N, sizeof(double));

    // Transpose A and B
    transpose(A, At, N);
    transpose(B, Bt, N);

    // Compute AtxB = A^T * B
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            register double aik = At[i * N + k];
            for (int j = 0; j < N; j++) {
                AtxB[i * N + j] += aik * B[k * N + j];
            }
        }
    }

    // Compute BxA = B * A
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            register double bik = B[i * N + k];
            for (int j = k; j < N; j++) {  // Only compute upper triangular part of A
                BxA[i * N + j] += bik * A[k * N + j];
            }
        }
    }

    // Sum AtxB and BxA
    for (int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }

    // Compute C = (AtxB + BxA) * B^T
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            register double tmp = AtxB_plus_BxA[i * N + k];
            for (int j = 0; j < N; j++) {
                C[i * N + j] += tmp * Bt[k * N + j];
            }
        }
    }

    free(At);
    free(Bt);
    free(AtxB);
    free(BxA);
    free(AtxB_plus_BxA);

    return C;
}