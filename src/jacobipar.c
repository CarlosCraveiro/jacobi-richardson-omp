#include <omp.h>
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define THRESHOLD 0.001

double max(const double* array, int size) {
    double greater = array[0];

    for(int i = 1; i < size; i++) {
        greater = (greater > array[i])? greater : array[i];    
    }

    return greater;
}

double gaussjacobi_error_parallel(const matrix_t* Xk, const matrix_t* Xkprev, int n_threads) {
    double *diff_vec = malloc(Xk->rows * sizeof(double));
    double *abs_xk_vec = malloc(Xk->rows * sizeof(double));
    
    #pragma omp parallel for num_threads(n_threads) shared(diff_vec, abs_xk_vec, Xk, Xkprev)
    for(int i = 0; i < Xk->rows; i++) {
        abs_xk_vec[i] = fabs(Xk->data[Xk->columns * i + 0]);
        diff_vec[i] = fabs(Xk->data[Xk->columns * i + 0] - Xkprev->data[Xkprev->columns * i + 0]);
    }
    
    double max_diff = max(diff_vec, Xk->rows); 
    double max_Xk_abs_element = max(abs_xk_vec, Xk->rows); 
    
    free(diff_vec);
    free(abs_xk_vec);

    return (max_diff / max_Xk_abs_element);
} 

matrix_t gaussjacobi_parallel(const matrix_t* A, const matrix_t* B, int n_threads) {
    matrix_t Xk = init_matrix(B->rows, 1, 1);
    matrix_t Xkprev = init_matrix(B->rows, 1, 1);
    
    #pragma omp parallel num_threads(n_threads) default(firstprivate) shared(A, B, Xk, Xkprev)
    { 
    #pragma omp single
    {
    //int itr = 0;
    do { 
        matrix_swap(&Xkprev, &Xk); /* FROM HERE */ 
        
        //#pragma omp parallel for num_threads(n_threads) default(private) shared(A, B, Xk, Xkprev)
        for(int i = 0; i < B->rows; i++) {
            #pragma omp task shared(A, B, Xk, Xkprev) firstprivate(i)
            {
            double xi = B->data[B->columns * i + 0];
            
            #pragma omp simd reduction(+:xi)
            for(int j = 0; j < A->columns ; j++) {
                xi += -1 * A->data[A->columns * i + j] * Xkprev.data[Xkprev.columns * j + 0];
            } 
            xi -= -1 * A->data[A->columns * i + i] * Xkprev.data[Xkprev.columns * i + 0];
            Xk.data[Xk.columns * i + 0] = xi / A->data[A->columns * i + i];
            }
        }
        #pragma omp taskwait
        //printf("Iteration %d\n", itr++);
        //print_matrix(&Xk, 1);
        /* TO HERE BECOME TASK 1*/
    } while(gaussjacobi_error_parallel(&Xk, &Xkprev, n_threads) > THRESHOLD);
    /* EVAL ERROR BECOME TASK NUMBER 2*/ // MAKE COPIES OF MATRIX
    }
    }
    free_matrix(Xkprev);
    return Xk;
}

int main(int argc, char* argv[]) {
    if(argc != 4) {
        printf("Incorrect number of arguments!\n");
        printf("Correct Usage:\n");
        printf("$ ./jacobipar <N> <T> <seed>\n");
        printf("\tN - Matrix order\n");
        printf("\tT - Number of threads\n");
        printf("\tseed - seed for the pseudorandom number generator\n");
        exit(-1);
    }
    int seed = atoi(argv[3]);
    int number_of_threads = atoi(argv[2]);
    int order = atoi(argv[1]);
    printf("test %d\n", seed);
    srand(seed);
    
    matrix_t A = init_rand_matrix(order, order);
    matrix_t B = init_rand_matrix(order, 1);
    
    //print_matrix(&A, 0);
    //print_matrix(&B, 0);

    matrix_t C = gaussjacobi_parallel(&A, &B, number_of_threads);
    
    printf("Result: \n");

    print_matrix(&C, 0);

    free_matrix(C);
    free_matrix(A);
    free_matrix(B);

    return 0;
 }
