#include <matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <math.h>

void matrix_swap(matrix_t* M1, matrix_t* M2) {
    long double* aux =  M1->data;
    M1->data = M2->data;
    M2->data = aux;
}

matrix_t matrix_transpose(const matrix_t* A) {
    matrix_t At;

    At.data = malloc(A->columns * A->rows * sizeof(long double)); // NxUninitialize 
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
            printf("%.3Lf", M->data[M->columns*i + j]);

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

matrix_t init_matrix(int rows, int columns, long double init_value) {
    // Definition of the matrix
    
    matrix_t A;

    A.data = malloc(rows * columns * sizeof(long double)); // NxUninitialize 
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

    A.data = malloc(rows * columns * sizeof(long double)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;
    long double diag = 0.0f;
    for(int i = 0; i < rows; i++) {  
        for(int j = 0; j < columns; j++) {
            long double gen_value = (rand() - (RAND_MAX/2)) / 1000.0f;
            A.data[A.columns * i + j] = gen_value;
            diag += fabsl(gen_value);
        }
        if(rows == columns) { // is a diagonal matrix
            A.data[A.columns * i + i] = diag - A.data[A.columns * i + i];
        } else if (columns == 1) {
            A.data[A.columns * 0 + i] *= rows * rows;
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
