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
 *        Código Ferramentas de Matrizes                     *
*************************************************************/

#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <stddef.h>
#include <math.h>

// Função para trocar os dados de duas matrizes
void matrix_swap(matrix_t* M1, matrix_t* M2) {
    matrix_value_t* aux =  M1->data;
    M1->data = M2->data;
    M2->data = aux;
}

// Função para transpor uma matriz
matrix_t matrix_transpose(const matrix_t* A) {
    matrix_t At;

    // Alocação de memória para a matriz transposta
    At.data = malloc(A->columns * A->rows * sizeof(*At.data)); // NxUninitialize 
    At.columns= A->columns;
    At.columns = A->rows; // Corrigido para At.rows

    // Transposição da matriz
    for(int i = 0; i < A->columns; i++) {
        for(int j = 0; j < A->rows; j++) {
            At.data[At.columns*i + j] = A->data[A->columns*i + j];
        }
    }

    return At;
}

// Função para inicializar uma matriz com um valor específico
matrix_t init_matrix(int rows, int columns, matrix_value_t init_value) {
    // Definição da matriz
    matrix_t A;

    // Alocação de memória para a matriz
    A.data = malloc(rows * columns * sizeof(*A.data)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;

    // Inicialização dos elementos da matriz
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < columns; j++) {
            A.data[i*A.columns + j] = init_value;
        }
    }

    return A;
}

// Função para inicializar uma matriz aleatória com dominância diagonal
matrix_t init_rand_diag_dominant_matrix(int order) {
    // Definição da matriz
    matrix_t A;

    // Alocação de memória para a matriz
    A.data = malloc(order * order * sizeof(*A.data)); // NxUninitialize 
    A.rows = order;
    A.columns = order;

    matrix_value_t diag = 0.0f;

    // Inicialização dos elementos da matriz
    for(int i = 0; i < order; i++) {  
        for(int j = 0; j < order; j++) {
            matrix_value_t gen_value = rand()%RANDOM_RANGE;
            A.data[A.columns * i + j] = gen_value;
            diag += gen_value;
        }

        // Garantir dominância diagonal
        A.data[A.columns * i + i] = diag - A.data[A.columns * i + i] + 1;
        diag = 0.0f;
    }
    
    return A;
}

// Função para inicializar uma matriz aleatória
matrix_t init_rand_matrix(int rows, int columns) {
    // Definição da matriz
    matrix_t A;

    // Alocação de memória para a matriz
    A.data = malloc(rows * columns * sizeof(*A.data)); // NxUninitialize 
    A.rows = rows;
    A.columns = columns;

    matrix_value_t diag = 0.0f;

    // Inicialização dos elementos da matriz
    for(int i = 0; i < rows; i++) {  
        for(int j = 0; j < columns; j++) {
            A.data[A.columns * i + j] = rand()%RANDOM_RANGE;
        }
    }
    
    return A;
}

// Função para calcular a i-ésima linha do produto matriz-vetor (GEMV)
matrix_value_t ith_row_GEMV(matrix_t* A, matrix_t* B, size_t row) {
    // Verifica se B é um vetor coluna
    if (B->columns != 1) {
        fprintf(stderr, "ith_row_GEMV: second argument must be columns vector\n");

    }
    
    // Verifica compatibilidade entre A e B
    if (A->columns != B->rows) {
        fprintf(stderr, "ith_row_GEMV: A and B must be compatible matrices\n");
    }

    matrix_value_t sum = 0;
    matrix_value_t error = 0;
    matrix_value_t temp = 0;
    matrix_value_t y = 0;

    // Aplicação da fórmula de soma compensada de Kahan para melhor precisão
    for(size_t j = 0; j < A->columns; j++) {
        temp = sum;
        y = (A->data)[A->columns*row + j] * (B->data)[j] + error;
        sum = temp + y;
        error = (temp - sum) + y;
    }

    return sum;
}

// Função para obter um elemento específico da matriz
matrix_value_t get_entry(matrix_t* A, size_t i, size_t j) {
    if (i >= A->rows || j >= A->columns) {
        fprintf(stderr, "get_entry: index pair (%zu, %zu) out of bounds.\n", i, j);
        return -1.0;
    }

    return (A->data)[A->columns*i + j];
}

// Função para liberar a memória de uma matriz
void free_matrix(matrix_t A) {
    free(A.data);
}