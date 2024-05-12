
#include <stdio.h>
#include <stdlib.h>

// C=(At×B+B×A)×Bt

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double *B) {
    printf("NEOPT SOLVER\n");
    
    double *C = calloc(N * N, sizeof(double));
    double *AtxB = calloc(N * N, sizeof(double));
    double *BxA = calloc(N * N, sizeof(double));

    // Calculate AtxB = A^T * B
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k <= i; k++) {  // A^T is lower triangular
                sum += A[k * N + i] * B[k * N + j];  // A^T[i][k] is A[k][i]
            }
            AtxB[i * N + j] = sum;
        }
    }

    printf ("AtxB\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%7.2f ", AtxB[i * N + j]);
        }
        printf("\n");
    }

    // Calculate BxA = B * A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k <= j; k++) {  
                sum += B[i * N + k] * A[k * N + j];
            }
            BxA[i * N + j] = sum;
        }
    }

    printf ("BxA\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%7.2f ", BxA[i * N + j]);
        }
        printf("\n");
    }

    // Calculate C = (AtxB + BxA) * Bt (final multiplication using B as B^T correctly)
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