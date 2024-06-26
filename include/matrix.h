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
 *       Header Código Ferramentas de Matrizes               *
*************************************************************/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stddef.h>

// Alcance dos números aleatórios gerados
#define RANDOM_RANGE 10

typedef double matrix_value_t;
typedef struct Matrix matrix_t;

// Estrutura de Dados de uma Matriz
struct Matrix {
    int rows;
    int columns;
    matrix_value_t* data;
};

// Declaração das funções

matrix_t matrix_transpose(const matrix_t* A);

void matrix_swap(matrix_t* M1, matrix_t* M2);

void matrix_copy(matrix_t* Copy, const matrix_t* Reference);

matrix_t init_rand_matrix(int rows, int columns);

matrix_t init_rand_diag_dominant_matrix(int order);

matrix_t init_matrix(int rows, int columns, matrix_value_t init_value); 

matrix_value_t ith_row_GEMV(matrix_t* A, matrix_t* B, size_t row);

matrix_value_t get_entry(matrix_t* A, size_t i, size_t j);

void free_matrix(matrix_t A);

#endif