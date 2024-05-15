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
 *             Código Jacobi Sequencial                      *
*************************************************************/

#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define THRESHOLD 0.001 // Define o limite de erro para a convergência

// Função para encontrar o maior valor em um array
matrix_value_t max(const matrix_value_t* array, int size) {
    matrix_value_t greater = array[0];

    for(int i = 1; i < size; i++) {
        greater = (greater > array[i])? greater : array[i];    
    }

    return greater;
}

// Função para calcular o erro do método de Jacobi
matrix_value_t gaussjacobi_error(const matrix_t* Xk, const matrix_t* Xkprev) {

    // Aloca vetores para armazenar as diferenças e os valores absolutos de Xk
    matrix_value_t *diff_vec = malloc(Xk->rows * sizeof(*diff_vec));
    matrix_value_t *abs_xk_vec = malloc(Xk->rows * sizeof(*abs_xk_vec));
    
    // Calcula a diferença absoluta entre Xk e Xkprev
    for(int i = 0; i < Xk->rows; i++) {
        abs_xk_vec[i] = fabs(Xk->data[Xk->columns * i + 0]);
        diff_vec[i] = fabs(Xk->data[Xk->columns * i + 0] - Xkprev->data[Xkprev->columns * i + 0]);
    }

    // Encontra os valores máximos dos vetores de diferença e valores absolutos
    matrix_value_t max_diff = max(diff_vec, Xk->rows); 
    matrix_value_t max_Xk_abs_element = max(abs_xk_vec, Xk->rows); 
    
    // Libera a memória dos vetores auxiliares
    free(diff_vec);
    free(abs_xk_vec);

    // Calcula e retorna o erro relativo
    return (max_diff / max_Xk_abs_element);
} 

// Função que implementa o método de Jacobi
matrix_t gaussjacobi(const matrix_t* A, const matrix_t* B) {

    // Inicializa os vetores Xk e Xkprev
    matrix_t Xk = init_matrix(B->rows, 1, 1);
    matrix_t Xkprev = init_matrix(B->rows, 1, 1);

    // Iteração do método de Jacobi
    do {
        // Troca os ponteiros de Xk e Xkprev
        matrix_swap(&Xkprev, &Xk); 
        
        // Atualiza cada elemento de Xk
        for(int i = 0; i < B->rows; i++) {
            matrix_value_t xi = B->data[B->columns * i + 0];

            // Soma ponderada dos elementos de A e Xkprev, excluindo a diagonal
            for(int j = 0; j < A->columns ; j++) {
                if(i != j) {
                    xi += -1 * A->data[A->columns * i + j] * Xkprev.data[Xkprev.columns * j + 0];
                }
            }

            // Atualiza o valor de Xk com o valor calculado
            Xk.data[Xk.columns * i + 0] = xi / A->data[A->columns * i + i];
        }

        // Continua até o erro ser menor que o limite
    } while(gaussjacobi_error(&Xk, &Xkprev) > THRESHOLD);
    
    // Libera a memória de Xkprev
    free_matrix(Xkprev);

    // Retorna a solução
    return Xk;
}

int main(int argc, char* argv[]) {
    // Verifica o número de argumentos
    if(argc != 4) {
        printf("Incorrect number of arguments!\n");
        printf("Correct Usage:\n");
        printf("$ ./jacobiseq <N> <seed> <row num>\n");
        printf("\tN - Matrix order\n");
        printf("\tseed - seed for the pseudorandom number generator\n");
        printf("\trow num - row in which the calculated B[i] is compared to its real value\n");
        exit(-1);
    }

    // Lê os argumentos de entrada
    int seed = atoi(argv[2]);
    int order = atoi(argv[1]);
    int row_index = atoi(argv[3]);

    srand(seed); // Seta a semente para o gerador de números aleatórios
    
    // Inicializa A com dominância diagonal e B aleatoriamente
    matrix_t A = init_rand_diag_dominant_matrix(order);
    matrix_t B = init_rand_matrix(order, 1);

    // Aplica o método de Jacobi
    matrix_t C = gaussjacobi(&A, &B);
    matrix_value_t calc_bi = ith_row_GEMV(&A, &C, (size_t)row_index);
    matrix_value_t real_bi = get_entry(&B, (size_t)row_index, 0);
    
    // Exibe os resultados
    printf("=============================================================\n");
    printf("Here are the results: \n");
    printf("\t- Real value of B[%zu] is: %f\n", row_index, real_bi);
    printf("\t- Calculated value of B[%zu] is: %f\n", row_index, calc_bi);
    printf("\t- Absolute error is: %f\n", fabs(real_bi - calc_bi));
    printf("\t- Relative error is: %f\n", fabs((real_bi - calc_bi)/real_bi));
    printf("=============================================================\n");

    // Libera a memória das matrizes
    free_matrix(C);
    free_matrix(A);
    free_matrix(B);

    return 0;
 }
