#include "utils.h"

// C = (At × B + B × A) × Bt
void transpose_matrices(double *A, double *transA, double *B, double *transB, int N)
{
	// Declare counters as registers
	register int i = 0, counter = 0;

	// Declare origin pointers
	double *At_origin = transA;
	double *Bt_origin = transB;

	// Declare iterators
	double *A_iterator = A;
	double *B_iterator = B;

	while (i++ < N)
	{
		// Start with iterators at the beginning of the matrices
		double *At_iterator = At_origin;
		double *Bt_iterator = Bt_origin;

		// Move origin pointers to the next row
		At_origin++;
		Bt_origin++;

		// Reset counters
		counter = 0;
		while (counter < N)
		{
			// Set the values of the transposed matrices
			*At_iterator = *A_iterator;
			*Bt_iterator = *B_iterator;

			// Move the transposed matrices iterators to the next element
			// +N to move to the next row
			At_iterator += N;
			Bt_iterator += N;

			// Move the iterators to the next element
			A_iterator++;
			B_iterator++;

			counter++;
		}
	}
}

double *fast_multiply(double *A, double *B, double *C, int N,
					  int start_i, int stop_i, int start_j, int stop_j,
					  int start_k, int stop_k, int triangular[4])
{
	// Declare origin pointers
	double *AB_origin = &C[0];
	double *B_origin = &B[0];
	double *A_origin = &A[0];

	// Calculate C = A * B
	for (register int i = start_i; i < stop_i; i++)
	{
		// A iterator starts at the beginning of the row
		register double *A_iterator = A_origin + i * N;

		// B iterator starts at the beginning of the matrix
		register double *B_iterator = B_origin;

		if (triangular[0] == 2)
		{ // If A is upper triangular
			A_iterator += i;
		}

		if (triangular[0] == 1)
		{ // If A is lower triangular
			stop_k = i + 1;
		}

		for (register int k = start_k; k < stop_k; k++)
		{
			// AB iterator starts at the beginning of the row
			register double *AB_iterator = AB_origin;
			// B interior iterator starts at the B current iterator
			register double *B_interior_iterator = B_iterator;
			// Just take the value of the current iterator for A
			register double A_interior_iterator = *A_iterator;

			if (triangular[1] == 2)
			{ // If B is upper triangular
				B_interior_iterator += k;
			}

			if (triangular[1] == 1)
			{ // If B is lower triangular
				stop_j = k + 1;
			}

			for (register int j = start_j; j < stop_j; j++)
			{
				// Place value in the AB matrix
				// AB[i][j] += A[i][k] * B[k][j]
				*AB_iterator += A_interior_iterator * *B_interior_iterator;

				// Move the result iterator to the next element
				AB_iterator++;

				// Move the B interior iterator to the next element
				B_interior_iterator++;
			}
			// Move the A iterator to the next element
			A_iterator++;

			// Move the B iterator to the next row
			B_iterator += N;
		}

		// Move the AB_origin to the next row
		AB_origin += N;
	}

	return C;
}

double *my_solver(int N, double *A, double *B)
{
	printf("OPT SOLVER\n");

	// Allocate memory for matrices
	double *C = calloc(N * N, sizeof(double));
	double *AtxB = calloc(N * N, sizeof(double));
	double *BxA = calloc(N * N, sizeof(double));
	double *At = calloc(N * N, sizeof(double));
	double *Bt = calloc(N * N, sizeof(double));
	double *AtxB_plus_BxA = calloc(N * N, sizeof(double));

	// Transpose A and B
	transpose_matrices(A, At, B, Bt, N);

	// Calculate At x B where At is lower triangular
	int triangular[2];
	triangular[0] = 1;
	triangular[1] = 0;
	AtxB = fast_multiply(At, B, AtxB, N, 0, N, 0, N, 0, N, triangular);

	// Calculate B * A where A is upper triangular
	triangular[0] = 0;
	triangular[1] = 0;
	BxA = fast_multiply(B, A, BxA, N, 0, N, 0, N, 0, N, triangular);

	// Calculate At x B + B x A
	for (int i = 0; i < N * N; i++)
	{
		AtxB_plus_BxA[i] = AtxB[i] + BxA[i];
	}

	// Calculate (At x B + B x A) * Bt
	triangular[0] = 0;
	triangular[1] = 0;
	C = fast_multiply(AtxB_plus_BxA, Bt, C, N, 0, N, 0, N, 0, N, triangular);

	// Free memory
	free(AtxB);
	free(BxA);
	free(At);
	free(Bt);
	free(AtxB_plus_BxA);

	return C;
}
