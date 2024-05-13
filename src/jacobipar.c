#include <omp.h>
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define THRESHOLD 0.001

matrix_value_t max(const matrix_value_t* array, int size) {
    matrix_value_t greater = array[0];

    for(int i = 1; i < size; i++) {
        greater = (greater > array[i])? greater : array[i];    
    }

    return greater;
}

matrix_value_t gaussjacobi_error_parallel(const matrix_t* Xk, const matrix_t* Xkprev) {
    matrix_value_t *diff_vec = malloc(Xk->rows * sizeof(*diff_vec));
    matrix_value_t *abs_xk_vec = malloc(Xk->rows * sizeof(*abs_xk_vec));
    
    #pragma omp simd
    for(int i = 0; i < Xk->rows; i++) {
        abs_xk_vec[i] = fabs(Xk->data[Xk->columns * i + 0]);
        diff_vec[i] = fabs(Xk->data[Xk->columns * i + 0] - Xkprev->data[Xkprev->columns * i + 0]);
    }
    
    matrix_value_t max_diff = max(diff_vec, Xk->rows); 
    matrix_value_t max_Xk_abs_element = max(abs_xk_vec, Xk->rows); 
    
    free(diff_vec);
    free(abs_xk_vec);

    return (max_diff / max_Xk_abs_element);
}

void gaussjacobi_iteration_parallel(const matrix_t* A, const matrix_t* B, const matrix_t* Xk, const matrix_t* Xkprev, int n_threads) {
    #pragma omp parallel for num_threads(n_threads)
    for(int i = 0; i < B->rows; i++) {
        matrix_value_t xi = 0.0;
        #pragma omp simd reduction(+:xi)
        for(int j = 0; j < A->columns; j++)
            xi += -1 * A->data[A->columns * i + j] * Xkprev->data[Xkprev->columns * j + 0];
         
	xi += B->data[B->columns * i + 0];
        xi += A->data[A->columns * i + i] * Xkprev->data[Xkprev->columns * i + 0]; //the diagonal is added cause of the previous loop
										 //so we remove it by adding witout multiplying with -1
        Xk->data[Xk->columns * i + 0] = xi / A->data[A->columns * i + i];
    }
}

matrix_t gaussjacobi_parallel(const matrix_t* A, const matrix_t* B, int n_threads) {
    matrix_t Xki = init_matrix(B->rows, 1, 1);    //current Xk (i) -- only 'iteration' will interact with this, so no lock
    matrix_t Xki_m1 = init_matrix(B->rows, 1, 1); //previous Xk (i-1) i_m1 -- will swap to the buffer of Xki at the end,
    matrix_t Xki_m2 = init_matrix(B->rows, 1, 1); //before the previous Xk (i-2)   i_m2 -- will swap to the buffer of Xki at the end,
    //the critical region only happens when it updates Xki_m1 and Xki_m2 by swapping the buffers at the end of 'iteration',
    //as they are treated as readonly on the other sections of the code, because of that,
    //only a single critical region called 'error-check-update' is needed, which is locked when 'iteration' updates 
    //the pointers at the end, and when 'error check' is reading both buffers
    //this is also used to signal 'error check' that 'iteration' ran, and that it can verify the stop condition with 'first_iter_ran'
    int continue_iter = 1;
    int iter_ran = 0;
    //1 thread for each task
    #pragma omp parallel num_threads(2) shared(A, B, Xki, Xki_m1, Xki_m2, continue_iter, iter_ran, n_threads)
    { 
        #pragma omp single nowait 
        {
            // The error check and the jacobi iterations will happen in different threads,
	    // where one can happen without waiting for the other
            #pragma omp task shared(Xki_m1, Xki_m2, continue_iter, iter_ran, n_threads) // gaussjacobi error check
            {
		int should_continue = 1;
                while(should_continue) {
                    #pragma omp critical (error_check_update)
            	    {
		        #pragma omp flush(iter_ran)
			if(iter_ran){ //only run if an iteration ran
            	            matrix_value_t error;
			    //half the threads for the parallel for
                            error = gaussjacobi_error_parallel(&Xki_m1, &Xki_m2); // Using the i-1 and i-2 to find the error
            	            if(error <= THRESHOLD){
            	                continue_iter = 0;
			        should_continue = 0;
			    }
			    iter_ran = 0; //marks that it read the last iteration
		            #pragma omp flush(iter_ran, continue_iter)
			}
            	    }
                    #pragma omp taskyield
                }
            }
            // gaussjacobi iteration
	    int should_continue = 1;
            while(should_continue) {
                gaussjacobi_iteration_parallel(A, B, &Xki, &Xki_m1, n_threads - 1);
                #pragma omp critical (error_check_update)
                {
                    matrix_swap(&Xki_m2, &Xki_m1); // Xki_m2 receives Xki_m1 buffer, and Xki_m1 receives "empty" buffer
                    matrix_swap(&Xki, &Xki_m1);    // Xki_m1 receives Xki buffer, and Xki receives "empty" buffer
            	    iter_ran = 1; // sets the flag that a iteration ran
		    #pragma omp flush(iter_ran, continue_iter)
		    if(continue_iter == 0){
			should_continue = 0;
		    }
            	}
            }
        }
    }
    
    free_matrix(Xki_m1);
    free_matrix(Xki_m2);
    return Xki;
}

int main(int argc, char* argv[]) { 
    omp_set_nested(1);
    omp_set_dynamic(0);

    if(argc != 5) {
        printf("Incorrect number of arguments!\n");
        printf("Correct Usage:\n");
        printf("$ ./jacobipar <N> <T> <seed> <row num>\n");
        printf("\tN - Matrix order\n");
        printf("\tT - Number of threads\n");
        printf("\tseed - seed for the pseudorandom number generator\n");
        printf("\trow num - row in which the calculated B[i] is compared to its real value\n");
        exit(-1);
    }

    int row_index = atoi(argv[4]);
    int seed = atoi(argv[3]);
    int number_of_threads = atoi(argv[2]);
    int order = atoi(argv[1]);
    printf("test %d\n", seed);
    srand(seed);
    
    matrix_t A = init_rand_diag_dominant_matrix(order);
    matrix_t B = init_rand_matrix(order, 1);
    
    //print_matrix(&A, 0);
    //print_matrix(&B, 0);

    matrix_t C = gaussjacobi_parallel(&A, &B, number_of_threads);
    matrix_value_t calc_bi = ith_row_GEMV(&A, &C, (size_t)row_index);
    matrix_value_t real_bi = get_entry(&B, (size_t)row_index, 0);
    
    printf("=============================================================\n");
    printf("Here are the results: \n");
    printf("\t- Real value of B[%zu] is: %f\n", row_index, real_bi);
    printf("\t- Calculated value of B[%zu] is: %f\n", row_index, calc_bi);
    printf("\t- Absolute error is: %f\n", fabs(real_bi - calc_bi));
    printf("\t- Relative error is: %f\n", fabs((real_bi - calc_bi)/real_bi));
    printf("=============================================================\n");

    free_matrix(C);
    free_matrix(A);
    free_matrix(B);

    return 0;
 }
