#include "utils.h"
#include <stdlib.h>
#include <stdio.h>

/*
 * Add your optimized implementation here
 */
void transpose_matrices(double *A, double *B, double *transA, double *transB, int N) {
    // Transpose A taking into account its upper triangular structure
    for (register int i = 0; i < N; i++) {
        double *ptr_A = &A[i * N];  // Pointer to the start of the current row of A
        double *ptr_transpose_A = &transA[i];  // Starting at the i-th position in transA
        
        for (register int j = 0; j <= i; j++) {
            *ptr_transpose_A = *ptr_A;  // Assign the element
            ptr_transpose_A += N;  // Move to the next row, same column in transA
            ptr_A++;  // Move to the next column in the current row of A
        }
    }

    // Transpose B normally
    for (register int i = 0; i < N; i++) {
        double *ptr_B = &B[i * N];  // Pointer to the start of the current row of B
        double *ptr_transpose_B = &transB[i];  // Starting at the i-th position in transB
        
        for (register int j = 0; j < N; j++) {
            *ptr_transpose_B = *ptr_B;  // Assign the element
            ptr_transpose_B += N;  // Move to the next row, same column in transB
            ptr_B++;  // Move to the next column in the current row of B
        }
    }
}

// Main solver function
double *my_solver(int N, double *A, double *B) {
    printf("OPT SOLVER\n");

    double *At = calloc(N * N, sizeof(double));
    double *Bt = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));
    double *AtxB_plus_BxA = calloc(N * N, sizeof(double));
    double *C = calloc(N * N, sizeof(double));

   	// Transpose A and B
    transpose_matrices(A, B, At, Bt, N);

    // Compute AtxB = A^T * B using more register hints
    for (register int i = 0; i < N; i++) {
        for (register int k = 0; k < N; k++) {
            register double aik = At[i * N + k];
            register double *AtxB_row = &AtxB[i * N];
            register double *B_row = &B[k * N];
            for (register int j = 0; j < N; j++) {
                AtxB_row[j] += aik * B_row[j];
            }
        }
    }

    // Compute BxA = B * A, leveraging A's upper triangular structure
    for (register int i = 0; i < N; i++) {
        for (register int k = 0; k < N; k++) {
            register double bik = B[i * N + k];
            register double *BxA_row = &BxA[i * N];
            register double *A_row = &A[k * N];
            for (register int j = k; j < N; j++) {
                BxA_row[j] += bik * A_row[j];
            }
        }
    }

    // Sum AtxB and BxA
    for (register int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }

    // Compute C = (AtxB + BxA) * B^T
    for (register int i = 0; i < N; i++) {
        register double *C_row = &C[i * N];
        for (register int k = 0; k < N; k++) {
            register double tmp = AtxB_plus_BxA[i * N + k];
            register double *Bt_row = &Bt[k * N];
            for (register int j = 0; j < N; j++) {
                C_row[j] += tmp * Bt_row[j];
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
