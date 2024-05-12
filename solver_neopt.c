/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"

// C=(At×B+B×A)×Bt

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	printf("NEOPT SOLVER\n");
	
	double *At = calloc(N * N, sizeof(double));
	double *AtxB = calloc(N * N, sizeof(double));
	double *BxA = calloc(N * N, sizeof(double));
	double *AtxB_plus_BxA = calloc(N * N, sizeof(double));
	double *C = calloc(N * N, sizeof(double));

	// Calculate At
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			At[i * N + j] = A[j * N + i];
		}
	}

	// Calculate AtxB
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				AtxB[i * N + j] += At[i * N + k] * B[k * N + j];
			}
		}
	}

	// Calculate BxA
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				BxA[i * N + j] += B[i * N + k] * A[k * N + j];
			}
		}
	}

	// Calculate AtxB_plus_BxA
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			AtxB_plus_BxA[i * N + j] = AtxB[i * N + j] + BxA[i * N + j];
		}
	}

	// Calculate C
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				C[i * N + j] += AtxB_plus_BxA[i * N + k] * B[j * N + k];
			}
		}
	}

	free(At);
	free(AtxB);
	free(BxA);
	free(AtxB_plus_BxA);

	return C;
}
