
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "cblas.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <string.h>

// C=(At×B+B×A)×Bt

void transpose_matrices(double *A, double *B, double *transA, double *transB,
						int N)
{
	double *ptr_orig_transpose_A = transA;
	double *ptr_orig_transpose_B = transB;

	double *ptr_A = A;
	double *ptr_B = B;
	register int i = 0;
	register int count;

	while (i < N)
	{
		double *ptr_transpose_A = ptr_orig_transpose_A;
		double *ptr_transpose_B = ptr_orig_transpose_B;

		count = 0;
		while (count < N)
		{
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

double *fast_multiply(double *A, double *B, double *C, int N,
					  int start_i, int stop_i,
					  int start_j, int stop_j,
					  int start_k, int stop_k, int triangular[4])
{
	double *pab_origin = &C[0];
	double *pb_origin = &B[0];
	double *pa_origin = &A[0];

	// Calculate C = A * B
	for (register int i = start_i; i < stop_i; i++)
	{
		register double *pa = pa_origin + i * N;
		register double *pb = pb_origin;

		if (triangular[0] == 2)
		{ // If A is upper triangular
			pa += i;
		}

		if (triangular[0] == 1)
		{ // If A is lower triangular
			stop_k = i + 1;
		}

		for (register int k = start_k; k < stop_k; k++)
		{
			register double *pab = pab_origin; // AB[i][j]
			register double *pb_k = pb;		   // &B[k * N]
			register double aik = *pa;

			if (triangular[1] == 2)
			{ // If A is upper triangular
				pb_k += k;
			}

			if (triangular[1] == 1)
			{ // If A is lower triangular
				stop_j = k + 1;
			}

			for (register int j = start_j; j < stop_j; j++)
			{
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

// BxA:
//    1.00    2.00   21.00
//    0.00   15.00   18.00
//    1.00    7.00   18.00

/*
 * Add your unoptimized implementation here
 */
// double* my_solver(int N, double *A, double *B) {
//     printf("OPT SOLVER\n");

//     double *C = calloc(N * N, sizeof(double));
//     double *AtxB = calloc(N * N, sizeof(double));
//     double *BxA = calloc(N * N, sizeof(double));
//     double *At = calloc(N * N, sizeof(double));
//     double *Bt = calloc(N * N, sizeof(double));
//     double *AtxB_plus_BxA = calloc(N * N, sizeof(double));
//     transpose_matrices(A, B, At, Bt, N);

// 	double *pab_origin = &AtxB[0];
// 	double *pb_origin = &B[0];
// 	double *pa_origin = &At[0];

// 	int triangular[2];
// 	triangular[0] = 1; // At is lower triangular
// 	triangular[1] = 0; // B is not triangular
// 	AtxB = fast_multiply(At, B, AtxB, N, 0, N, 0, N, 0, N, triangular);

//     // print the intermediate matrix AtxB
//     printf("AtxB:\n");
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             printf("%7.2f ", AtxB[i * N + j]);
//         }
//         printf("\n");
//     }

// 	triangular[0] = 0; // B is not triangular
// 	triangular[1] = 0; // A is upper triangular
// 	BxA = fast_multiply(B, A, BxA, N, 0, N, 0, N, 0, N, triangular);

//     // Calculate AtxB + BxA
//     for (int i = 0; i < N * N; i++) {
//         AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
//     }

//     // PRINT BxA
//     printf("BxA:\n");
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             printf("%7.2f ", BxA[i * N + j]);
//         }
//         printf("\n");
//     }

// 	triangular[0] = 0;
// 	triangular[1] = 0;
// 	C = fast_multiply(AtxB_plus_BxA, Bt, C, N, 0, N, 0, N, 0, N, triangular);

//     free(AtxB);
//     free(BxA);

//     return C;
// }

double *my_solver(int N, double *A, double *B)
{
	printf("BLAS SOLVER\n");

	double *C = calloc(N * N, sizeof(double));
	double *B_copy = calloc(N * N, sizeof(double));
	double *B_copy1 = calloc(N * N, sizeof(double));

	// Copy B to B_copy
	cblas_dcopy(N * N, B, 1, B_copy, 1);
	cblas_dcopy(N * N, B, 1, B_copy1, 1);

	// Calculate A^T * B where A^T is lower triangular
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit, N, N, 1.0, A, N, B, N);

	// Print AtxB
	printf("AtxB:\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%7.2f ", B[i * N + j]);
		}
		printf("\n");
	}

	// Calculate B_copy * A where A is upper triangular
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1.0, A, N, B_copy, N);

	// Print BxA
	printf("BxA:\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%7.2f ", B_copy[i * N + j]);
		}
		printf("\n");
	}

	// Sum AtxB and BxA using cblas_daxpy
	// Where AtxB is stored in B and BxA is stored in B_copy
	cblas_daxpy(N * N, 1.0, B_copy, 1, B, 1);

	// Calculate (AtxB + BxA) * B^T using BxA which now holds the sum of AtxB and BxA
	// Where B^T is stored in B_copy1 and (AtxB + BxA) is stored in B
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, B, N, B_copy1, N, 0.0, C, N);

	free(B_copy);
	free(B_copy1);

	return C;
}

// Function to print matrices
void print_matrix(const char *name, double *matrix, int N)
{
	printf("%s:\n", name);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%7.2f ", matrix[i * N + j]);
		}
		printf("\n");
	}
}

// Testing function
void test_my_solver()
{
	int N = 3;	  // Matrix dimension
	double A[] = {// Upper triangular matrix
				  1, 2, 3,
				  0, 5, 6,
				  0, 0, 9};
	double B[] = {// General matrix
				  1, 0, 2,
				  0, 3, 0,
				  1, 1, 1};

	// Allocate memory for the resulting matrix C
	double *C = my_solver(N, A, B);

	// Print the input matrices and the result
	print_matrix("Matrix A", A, N);
	print_matrix("Matrix B", B, N);
	print_matrix("Matrix C (Result)", C, N);

	// Free allocated memory for matrix C
	free(C);
}

int main()
{
	test_my_solver();
	return 0;
}