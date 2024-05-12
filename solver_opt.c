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

double* my_solver(int N, double *A, double *B) {
    printf("OPT SOLVER\n");
    
    double *C = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));

	double *pab_origin = &AtxB[0];
	double *pb_origin = &B[0];
	double *pa_origin = &A[0];

    // Calculate AtxB = A^T * B
	for (register int i = 0; i < N; i++) {
		register double *pb_k = pb_origin; // &B[i * N]
		register double *pa_k =
			pa_origin + i; // Sup triangular matrix, start from i
		for (register int k = i; k < N; k++) {
			register double *pab_j = pab_origin; // A[i][0]
			register double *pb_j = pb_k;		 // B[k][0]
			register double aik = *pa_k;		 // A[i * N + k];
			for (register int j = 0; j < N; j++) {
				*pab_j += aik * *pb_j;
				pab_j++; // AB[i][j++];
				pb_j++;	 // B[k][j++]
			}
			pb_k += N; // B[k++][j]
			pa_k++;
		}
		pab_origin += N; // AB[i++][j]
		pa_origin += N;	 // A[i++][k]
		pb_origin += N;	 // reset origin to next line
	}


    // Calculate BxA = B * A
    for (register int i = 0; i < N; i++) {
        register double sum = 0.0;
		register double *ptr_A = A;
		register double *ptr_B = B + i;
		for (register int j = 0; j < N; j++) {
			register double *ptr_A_k = ptr_A;
			register double *ptr_B_k = ptr_B;
			for (register int k = 0; k <= j; k++) {
				sum += *ptr_B_k * *ptr_A_k;
				ptr_A_k++;
				ptr_B_k++;
			}
			BxA[i * N + j] = sum;
			sum = 0.0;
			ptr_A += N;
		}
	}

    // Calculate C = (AtxB + BxA) * Bt (final multiplication using B as B^T correctly)
    for (register int i = 0; i < N; i++) {
        register double *ptr_B = B;
		for (register int j = 0; j < N; j++) {
			register double sum = 0.0;
			for (register int k = 0; k < N; k++) {
				sum += (AtxB[i * N + k] + BxA[i * N + k]) * *ptr_B;
				ptr_B++;
			}
			C[i * N + j] = sum;
		}
	}

    free(AtxB);
    free(BxA);

    return C;
}