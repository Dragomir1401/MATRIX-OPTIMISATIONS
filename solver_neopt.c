/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"

// C=(At×B+B×A)×Bt

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double *B) {
    printf("NEOPT SOLVER\n");
    
    double *C = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));

    // Calculate At x B
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k <= i; k++) {  // A^T is lower triangular
                sum += A[k * N + i] * B[k * N + j];  // A^T[i][k] is A[k][i]
            }
            AtxB[i * N + j] = sum;
        }
    }

    // Calculate B x A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k <= j; k++) {  
                sum += B[i * N + k] * A[k * N + j];
            }
            BxA[i * N + j] = sum;
        }
    }

    // Calculate C = (AtxB + BxA) x Bt (final multiplication using B as B^T correctly)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C[i * N + j] += (AtxB[i * N + k] + BxA[i * N + k]) * B[j * N + k];
            }
        }
    }

    free(AtxB);
    free(BxA);

    return C;
}