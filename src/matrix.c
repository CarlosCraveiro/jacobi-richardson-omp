#include <matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>

// Deprecated!!
void matrix_copy(matrix_t* Copy, const matrix_t* Reference) {
    for(int i = 0; i < Copy->lines; i++) {
        for(int j = 0; j < Copy->columns; j++) {
            Copy->data[i][j] = Reference->data[i][j];
        }
    }
}

void matrix_swap(matrix_t* M1, matrix_t* M2) {
    double** aux =  M1->data;
    M1->data = M2->data;
    M2->data = aux;
}

void print_matrix(matrix_t* M, int is_array) {
     
    if(is_array == 1) {
         printf("[ ");
    }

    for(int i = 0; i < M->lines; i++) {
        for(int j = 0; j < M->columns; j++) {
            printf("%.3lf", M->data[i][j]);

            if(j != (M->columns - 1) ||
                    ((is_array == 1) && (i != M->lines - 1))
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

matrix_t init_matrix(int lines, int columns, double init_value) {
    // Definition of the matrix
    
    matrix_t A;

    A.data = malloc(lines * sizeof(double *)); // NxUninitialize 
    A.lines = lines;
    A.columns = columns;
    for(int i = 0; i < lines; i++) {
        A.data[i] = malloc(lines * sizeof(double)); // NxN
        for(int j = 0; j < columns; j++) {
            A.data[i][j] = init_value;
        }
    }
    return A;
}

matrix_t init_rand_matrix(int lines, int columns) {
    // Definition of the matrix
    
    matrix_t A;

    A.data = malloc(lines * sizeof(double *)); // NxUninitialize 
    A.lines = lines;
    A.columns = columns;
    for(int i = 0; i < lines; i++) {
        A.data[i] = malloc(lines * sizeof(double)); // NxN
        int diag_int = rand();
        int upper_limit = diag_int / MAX(lines, columns);
        for(int j = 0; j < columns; j++) {
            if(i == j) {
                A.data[i][j] = (diag_int / 1000.0f);
            } else {
                double gen_value = ((rand() % (2*upper_limit)) - upper_limit) / 1000.0f;
                A.data[i][j] = gen_value;
            } 
        }
    }
    return A;
}

matrix_t multiply_matrices(matrix_t* A, matrix_t* B) {
    // Defines result matrix C
    matrix_t C = init_matrix(A->lines, B->columns, 0.0);
    
    for(int i = 0; i < C.lines; i++) {
        for(int j = 0; j < C.columns; j++) {
            for(int k = 0; k < C.lines; k++) { // Scalar product
                C.data[i][j] += A->data[i][k] * B->data[k][j];
            }
        }
    }

    return C;
}

void free_matrix(matrix_t A) {
    for (int i = 0; i < A.lines; i++) {
        free(A.data[i]);
    }
    free(A.data);
}
