#include "utils.h"
#include <stdlib.h>
#include <stdio.h>

// C = (At × B + B × A) × Bt

/*
 * Add your optimized implementation here
 */
// Helper function to transpose matrices
void transpose_matrices(double *A, double *B, double *transA, double *transB,
					  int N) {
	double *ptr_orig_transpose_A = transA;
	double *ptr_orig_transpose_B = transB;

	double *ptr_A = A;
	double *ptr_B = B;
	register int i = 0;
	register int count;

	while (i < N) {
		double *ptr_transpose_A = ptr_orig_transpose_A;
		double *ptr_transpose_B = ptr_orig_transpose_B;

		count = 0;
		while (count < N) {
			*ptr_transpose_A = *ptr_A;
			ptr_transpose_A += N;
			ptr_A++;

			*ptr_transpose_B = *ptr_B;
			ptr_transpose_B += N;
			ptr_B++;

			count++;
		}
		i++;
		ptr_orig_transpose_A++;
		ptr_orig_transpose_B++;
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
    transpose_matrices(A, B, At, Bt, N);

    // Compute AtxB = At * B
    for (register int i = 0; i < N; i++) {
        double register *At_row = At + i * N;
        double register *AtxB_row = AtxB + i * N;
        for (register int j = 0; j < N; j++) {
            double register sum = 0.0;
            double register *B_col = B + j;
            for (register int k = 0; k < N; k++) {
                sum += At_row[k] * B_col[k * N];
            }
            AtxB_row[j] = sum;
        }
    }

    // Compute BxA = B * A
    for (register int i = 0; i < N; i++) {
        register double *B_row = B + i * N;
        register double *BxA_row = BxA + i * N;
        for (register int j = 0; j < N; j++) {
            register double sum = 0.0;
            register double *A_col = A + j;
            for (register int k = 0; k < N; k++) {
                sum += B_row[k] * A_col[k * N];
            }
            BxA_row[j] = sum;
        }
    }

    // Compute AtxB_plus_BxA = AtxB + BxA
    for (int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }

    // Compute C = (AtxB_plus_BxA) * Bt
    for (register int i = 0; i < N; i++) {
        register double *AtxB_plus_BxA_row = AtxB_plus_BxA + i * N;
        register double *C_row = C + i * N;
        for (register int j = 0; j < N; j++) {
            register double sum = 0.0;
            register double *Bt_col = Bt + j;
            for (register int k = 0; k < N; k++) {
                sum += AtxB_plus_BxA_row[k] * Bt_col[k * N];
            }
            C_row[j] = sum;
        }
    }

    free(At);
    free(Bt);
    free(AtxB);
    free(BxA);
    free(AtxB_plus_BxA);

    return C;
}