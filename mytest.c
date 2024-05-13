
#include <stdio.h>
#include <stdlib.h>

// C=(At×B+B×A)×Bt



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

	double *pab_origin = &AtxB[0];
	double *pb_origin = &B[0];
	double *pa_origin = &At[0];

    // // Calculate AtxB = A^T * B
    // for (register int i = 0; i < N; i++) {
	// 	register double *pa = pa_origin + i * N;
	// 	register double *pb = pb_origin;
	// 	for (register int k = 0; k <= i; k++) { // A^T is lower triangular
	// 		register double *pab = pab_origin; // AB[i][j]
	// 		register double *pb_k = pb;		   // &B[k * N]
	// 		register double aik = *pa;
	// 		for (register int j = 0; j < N; j++) {
	// 			*pab += aik * *pb_k; // AB[i][j] = A[i][k] + B[k][j]
	// 			pab++;				 // AB[i][j++]
	// 			pb_k++;				 // B[k][j++]
	// 		}
	// 		pa++;	 // A[i][k++]
	// 		pb += N; // B[k++][j]
	// 	}
	// 	pab_origin += N;
	// }

	AtxB = fast_multiply(At, B, AtxB, N, 0, N, 0, N, 0, N);

    // print the intermediate matrix AtxB
    printf("AtxB:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%7.2f ", AtxB[i * N + j]);
        }
        printf("\n");
    }


	// pab_origin = &BxA[0];
	// pb_origin = &A[0];
	// pa_origin = &B[0];

    // // Calculate B x A = B * A
    // for (register int i = 0; i < N; i++) {
	// 	register double *pa = pa_origin + i * N;
	// 	register double *pb = pb_origin;
	// 	for (register int k = 0; k < N; k++) {
	// 		register double *pab = pab_origin; // AB[i][j]
	// 		register double *pb_k = pb;		   // &A[k * N]
	// 		register double aik = *pa;
	// 		for (register int j = 0; j < N; j++) {
	// 			*pab += aik * *pb_k; // AB[i][j] = A[i][k] + B[k][j]
	// 			pab++;				 // AB[i][j++]
	// 			pb_k++;				 // B[k][j++]
	// 		}
	// 		pa++;	 // A[i][k++]
	// 		pb += N; // B[k++][j]
	// 	}
	// 	pab_origin += N;
	// }

	BxA = fast_multiply(B, A, BxA, N, 0, N, 0, N, 0, N);

// 	BxA:
//    1.00    2.00   21.00 
//    0.00   15.00   18.00 
//    1.00    7.00   18.00 

    // Calculate AtxB + BxA
    for (int i = 0; i < N * N; i++) {
        AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
    }


    // PRINT BxA
    printf("BxA:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%7.2f ", BxA[i * N + j]);
        }
        printf("\n");
    }

    // pab_origin = &C[0];
	// pb_origin = &Bt[0];
	// pa_origin = &AtxB_plus_BxA[0];

    // // Calculate C = (AtxB + BxA) * Bt (final multiplication using B as B^T correctly)
    // for (register int i = 0; i < N; i++) {
	// 	register double *pa = pa_origin + i * N;
	// 	register double *pb = pb_origin;
	// 	for (register int k = 0; k < N; k++) {
	// 		register double *pab = pab_origin; // AB[i][j]
	// 		register double *pb_k = pb;		   // &B[k * N]
	// 		register double aik = *pa;
	// 		for (register int j = 0; j < N; j++) {
	// 			*pab += aik * *pb_k; // AB[i][j] = A[i][k] + B[k][j]
	// 			pab++;				 // AB[i][j++]
	// 			pb_k++;				 // B[k][j++]
	// 		}
	// 		pa++;	 // A[i][k++]
	// 		pb += N; // B[k++][j]
	// 	}
	// 	pab_origin += N;
	// }

	C = fast_multiply(AtxB_plus_BxA, Bt, C, N, 0, N, 0, N, 0, N);

    free(AtxB);
    free(BxA);

    return C;
}

// Function to print matrices
void print_matrix(const char* name, double *matrix, int N) {
    printf("%s:\n", name);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%7.2f ", matrix[i * N + j]);
        }
        printf("\n");
    }
}

// Testing function
void test_my_solver() {
    int N = 3;  // Matrix dimension
    double A[] = {  // Upper triangular matrix
        1, 2, 3,
        0, 5, 6,
        0, 0, 9
    };
    double B[] = {  // General matrix
        1, 0, 2,
        0, 3, 0,
        1, 1, 1
    };

    // Allocate memory for the resulting matrix C
    double *C = my_solver(N, A, B);
    
    // Print the input matrices and the result
    print_matrix("Matrix A", A, N);
    print_matrix("Matrix B", B, N);
    print_matrix("Matrix C (Result)", C, N);

    // Free allocated memory for matrix C
    free(C);
}

int main() {
    test_my_solver();
    return 0;
}