#include <matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>

void matrix_swap(matrix_t* M1, matrix_t* M2) {
    double* aux =  M1->data;
    M1->data = M2->data;
    M2->data = aux;
}

matrix_t matrix_transpose(const matrix_t* A) {
    matrix_t At;

    At.data = malloc(A->columns * A->rows * sizeof(double)); // NxUninitialize 
    At.columns= A->columns;
    At.columns = A->rows;
    for(int i = 0; i < A->columns; i++) {
        for(int j = 0; j < A->rows; j++) {
            At.data[At.columns*i + j] = A->data[A->columns*i + j];
        }
    }

    return At;
}

void print_matrix(matrix_t* M, int is_array) {
     
    if(is_array == 1) {
         printf("[ ");
    }

    for(int i = 0; i < M->rows; i++) {
        for(int j = 0; j < M->columns; j++) {
            printf("%.3lf", M->data[M->columns*i + j]);

            if(j != (M->columns - 1) ||
                    ((is_array == 1) && (i != M->columns - 1))
            ) {
                printf(", ");
            }
        }
        if(is_array == 0) {
            printf("\n");
        }
    }

    if(is_array == 1) {
        printf(" ]\n");
    } 
}

matrix_t init_matrix(int rows, int columns, double init_value) {
    // Definition of the matrix
    
    matrix_t A;

    A.data = malloc(rows * columns * sizeof(double)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < columns; j++) {
            A.data[i*A.columns + j] = init_value;
        }
    }
    return A;
}

matrix_t init_rand_matrix(int rows, int columns) {
    // Definition of the matrix
    
    matrix_t A;

    A.data = malloc(rows * columns * sizeof(double)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;
    for(int i = 0; i < rows; i++) {
        int diag_int = rand();
        int upper_limit = diag_int / MAX(rows, columns);
        for(int j = 0; j < columns; j++) {
            if(i == j) {
                A.data[A.columns*i + j] = (diag_int / 1000.0f);
            } else {
                double gen_value = ((rand() % (2*upper_limit)) - upper_limit) / 1000.0f;
                A.data[A.columns * i + j] = gen_value;
            } 
        }
    }
    return A;
}

matrix_t multiply_matrices(matrix_t* A, matrix_t* B) {
    // Defines result matrix C
    matrix_t C = init_matrix(A->rows, B->columns, 0.0);
    
    for(int i = 0; i < C.rows; i++) {
        for(int j = 0; j < C.columns; j++) {
            for(int k = 0; k < C.rows; k++) { // Scalar product
                C.data[i*C.columns + j] += A->data[i*A->columns + k] * B->data[B->columns * k + j];
            }
        }
    }

    return C;
}

void free_matrix(matrix_t A) {
    free(A.data);
}
