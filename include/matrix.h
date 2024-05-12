#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stddef.h>

#define RANDOM_RANGE 10

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

matrix_t init_rand_diag_dominant_matrix(int order);

matrix_t init_matrix(int rows, int columns, matrix_value_t init_value); 

matrix_t multiply_matrices(matrix_t* A, matrix_t* B);

matrix_value_t ith_row_GEMV(matrix_t* A, matrix_t* B, size_t row);

matrix_value_t get_entry(matrix_t* A, size_t i, size_t j);

void free_matrix(matrix_t A);

#endif