/************************************************************
 *          Trabalho Prático 1 - SCC0903                     *
 *                                                           *
 *      Nome: Artur Brenner Weber                            *
 *      nUSP: 12675451                                       *
 *      Nome: Carlos Henrique Craveiro Aquino Veras          *
 *      nUSP: 12547187                                       *
 *      Nome: Gabriel Franceschi Libardi                     *
 *      nUSP: 11760739                                       *
 *      Nome: Ivan Roberto Pancheniak                        *
 *      nUSP: 12624224                                       *
 *      Data de última atualizacao: 15/5/2024                *
 *      Ambiente: VSCode 1.89.1                              *
 *                                                           *
*************************************************************/

#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <stddef.h>
#include <math.h>

void matrix_swap(matrix_t* M1, matrix_t* M2) {
    matrix_value_t* aux =  M1->data;
    M1->data = M2->data;
    M2->data = aux;
}

matrix_t matrix_transpose(const matrix_t* A) {
    matrix_t At;

    At.data = malloc(A->columns * A->rows * sizeof(*At.data)); // NxUninitialize 
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
            printf("%.3f", M->data[M->columns*i + j]);

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

matrix_t init_matrix(int rows, int columns, matrix_value_t init_value) {
    // Definition of the matrix
    matrix_t A;

    A.data = malloc(rows * columns * sizeof(*A.data)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < columns; j++) {
            A.data[i*A.columns + j] = init_value;
        }
    }

    return A;
}

matrix_t init_rand_diag_dominant_matrix(int order) {
    // Definition of the matrix
    matrix_t A;

    A.data = malloc(order * order * sizeof(*A.data)); // NxUninitialize 
    A.rows = order;
    A.columns = order;

    matrix_value_t diag = 0.0f;
    for(int i = 0; i < order; i++) {  
        for(int j = 0; j < order; j++) {
            matrix_value_t gen_value = rand()%RANDOM_RANGE;
            A.data[A.columns * i + j] = gen_value;
            diag += gen_value;
        }

        A.data[A.columns * i + i] = diag - A.data[A.columns * i + i] + 1;
        diag = 0.0f;
    }
    
    return A;
}

matrix_t init_rand_matrix(int rows, int columns) {
    // Definition of the matrix
    matrix_t A;

    A.data = malloc(rows * columns * sizeof(*A.data)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;

    matrix_value_t diag = 0.0f;
    for(int i = 0; i < rows; i++) {  
        for(int j = 0; j < columns; j++) {
            A.data[A.columns * i + j] = rand()%RANDOM_RANGE;
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

matrix_value_t ith_row_GEMV(matrix_t* A, matrix_t* B, size_t row) {
    if (B->columns != 1) {
        fprintf(stderr, "ith_row_GEMV: second argument must be columns vector\n");

    }
    
    if (A->columns != B->rows) {
        fprintf(stderr, "ith_row_GEMV: A and B must be compatible matrices\n");
    }

    matrix_value_t sum = 0;
    matrix_value_t error = 0;
    matrix_value_t temp = 0;
    matrix_value_t y = 0;

    // Apply Kahan's compensated summation formula for better accuracy.
    for(size_t j = 0; j < A->columns; j++) {
        temp = sum;
        y = (A->data)[A->columns*row + j] * (B->data)[j] + error;
        sum = temp + y;
        error = (temp - sum) + y;
    }

    return sum;
}

matrix_value_t get_entry(matrix_t* A, size_t i, size_t j) {
    if (i >= A->rows || j >= A->columns) {
        fprintf(stderr, "get_entry: index pair (%zu, %zu) out of bounds.\n", i, j);
        return -1.0;
    }

    return (A->data)[A->columns*i + j];
}

void free_matrix(matrix_t A) {
    free(A.data);
}
