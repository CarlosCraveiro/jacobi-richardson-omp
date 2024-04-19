#pragma once

typedef struct Matrix matrix_t;

struct Matrix {
    int lines;
    int columns;
    double** data;
};

void matrix_copy(matrix_t* Copy, const matrix_t* Reference);

void print_matrix(matrix_t* M, int is_array);

matrix_t init_matrix(int lines, int columns, double init_value); 

matrix_t multiply_matrices(matrix_t* A, matrix_t* B);

void free_matrix(matrix_t A);
