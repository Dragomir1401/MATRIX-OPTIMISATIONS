#include "utils.h"

// C = (At × B + B × A) × Bt
void transpose_matrices(double *A, double *B, double *transA, double *transB, int N) {
	register int i = 0, counter = 0;
	double *At_origin = transA;
	double *A_iterator = A;
	double *Bt_origin = transB;
	double *B_iterator = B;

	while (i < N) {
		double *At_iterator = At_origin;
		double *Bt_iterator = Bt_origin;
		At_origin++;
		Bt_origin++;

		counter = 0;
		while (counter < N) {
			*At_iterator = *A_iterator;
			At_iterator += N;
			A_iterator++;

			*Bt_iterator = *B_iterator;
			Bt_iterator += N;
			B_iterator++;

			counter++;
		}

		i++;
	}
}

double* fast_multiply(double *A, double *B, double* C, int N,
					int start_i, int stop_i, 
					int start_j, int stop_j,
					int start_k, int stop_k, int triangular[4]) {
	double *pab_origin = &C[0];
	double *pb_origin = &B[0];
	double *pa_origin = &A[0];

    // Calculate C = A * B
    for (register int i = start_i; i < stop_i; i++) {
		register double *pa = pa_origin + i * N;
		register double *pb = pb_origin;

		if (triangular[0] == 2) { // If A is upper triangular
			pa += i;
		}

		if (triangular[0] == 1) { // If A is lower triangular
			stop_k = i + 1;
		}

		for (register int k = start_k; k < stop_k; k++) {
			register double *pab = pab_origin; // AB[i][j]
			register double *pb_k = pb;		   // &B[k * N]
			register double aik = *pa;

			if (triangular[1] == 2) { // If B is upper triangular
				pb_k += k;
			}

			if (triangular[1] == 1) { // If B is lower triangular
				stop_j = k + 1;
			}

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

double* my_solver(int N, double *A, double *B) {
    printf("OPT SOLVER\n");
    
    double *C = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));
    double *At = calloc(N * N, sizeof(double));
    double *Bt = calloc(N * N, sizeof(double));
    double *AtxB_plus_BxA = calloc(N * N, sizeof(double));
    transpose_matrices(A, B, At, Bt, N);

	int triangular[2];
	triangular[0] = 1;
	triangular[1] = 0;
	AtxB = fast_multiply(At, B, AtxB, N, 0, N, 0, N, 0, N, triangular);


	triangular[0] = 0;
	triangular[1] = 0;
	BxA = fast_multiply(B, A, BxA, N, 0, N, 0, N, 0, N, triangular);

    // Calculate AtxB + BxA
    for (int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }

	triangular[0] = 0;
	triangular[1] = 0;
	C = fast_multiply(AtxB_plus_BxA, Bt, C, N, 0, N, 0, N, 0, N, triangular);

    free(AtxB);
    free(BxA);
	free(At);
	free(Bt);
	free(AtxB_plus_BxA);

    return C;
}