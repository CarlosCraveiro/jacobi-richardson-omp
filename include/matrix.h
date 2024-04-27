#pragma once

typedef double matrix_value_t;
typedef struct Matrix matrix_t;

struct Matrix {
    int rows;
    int columns;
    matrix_value_t* data;
};

matrix_t matrix_transpose(const matrix_t* A);

void matrix_swap(matrix_t* M1, matrix_t* M2);

void matrix_copy(matrix_t* Copy, const matrix_t* Reference);

matrix_t init_rand_matrix(int rows, int columns);

void print_matrix(matrix_t* M, int is_array);

matrix_t init_matrix(int rows, int columns, matrix_value_t init_value); 

matrix_t multiply_matrices(matrix_t* A, matrix_t* B);

void free_matrix(matrix_t A);
