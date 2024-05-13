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

double* fast_multiply(double *A, double *B, double* C, int N,
					int start_i, int stop_i, 
					int start_j, int stop_j,
					int start_k, int stop_k) {
	double *pab_origin = &C[0];
	double *pb_origin = &B[0];
	double *pa_origin = &A[0];

    // Calculate C = A * B
    for (register int i = start_i; i < stop_i; i++) {
		register double *pa = pa_origin + i * N;
		register double *pb = pb_origin;

		for (register int k = start_k; k < stop_k; k++) {
			register double *pab = pab_origin; // AB[i][j]
			register double *pb_k = pb;		   // &B[k * N]
			register double aik = *pa;

			for (register int j = start_j; j < stop_j; j++) {
				*pab += aik * *pb_k; // AB[i][j] = A[i][k] + B[k][j]
				pab++;				 // AB[i][j++]
				pb_k++;				 // B[k][j++]
			}
			pa++;	 // A[i][k++]
			pb += N; // B[k++][j]
		}
		pab_origin += N;
	}

	return C;
}

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double *B) {
    printf("OPT SOLVER\n");
    
    double *C = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));
    double *At = calloc(N * N, sizeof(double));
    double *Bt = calloc(N * N, sizeof(double));
    double *AtxB_plus_BxA = calloc(N * N, sizeof(double));
    transpose_matrices(A, B, At, Bt, N);

	AtxB = fast_multiply(At, B, AtxB, N, 0, N, 0, N, 0, N);

	BxA = fast_multiply(B, A, BxA, N, 0, N, 0, N, 0, N);

    // Calculate AtxB + BxA
    for (int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }

	C = fast_multiply(AtxB_plus_BxA, Bt, C, N, 0, N, 0, N, 0, N);

    free(AtxB);
    free(BxA);
	free(At);
	free(Bt);
	free(AtxB_plus_BxA);

    return C;
}
