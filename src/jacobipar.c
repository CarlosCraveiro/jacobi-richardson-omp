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

matrix_t* gaussjacobi_parallel(const matrix_t* A, const matrix_t* B, int n_threads) {
    matrix_t *Xki = init_matrix(B->rows, 1, 1); //current Xk (i) -- only 'iteration' will interact with this, so no lock
    matrix_t *Xki_m1 = init_matrix(B->rows, 1, 1); //previous Xk (i-1) i_m1 -- will swap to the buffer of Xki at the end,
    matrix_t *Xki_m2 = init_matrix(B->rows, 1, 1); //before the previous Xk (i-2)   i_m2 -- will swap to the buffer of Xki at the end,
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
                            error = gaussjacobi_error_parallel(Xki_m1, Xki_m2); // Using the i-1 and i-2 to find the error
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
                gaussjacobi_iteration_parallel(A, B, Xki, Xki_m1, n_threads - 1);
                #pragma omp critical (error_check_update)
                {
                    matrix_swap(Xki_m2, Xki_m1); // Xki_m2 receives Xki_m1 buffer, and Xki_m1 receives "empty" buffer
                    matrix_swap(Xki, Xki_m1);    // Xki_m1 receives Xki buffer, and Xki receives "empty" buffer
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

    int order, seed, number_of_threads;
    int quiet_mode = 0;
    int incorrect_usage_of_args = 0;
    int number_argument = 0;
    for(int i = 1; i < argc; i++) {
        if(argv[i][0] == '-') { // this is a flag
	    if(strcmp(argv[i],"-q") == 0) quiet_mode = 1; //quiet flag
	    else if(strcmp(argv[i],"-h") == 0) incorrect_usage_of_args = 1; //force help text
	    else if(strcmp(argv[i],"-v") == 0) printf("Version: CAD\n");
            else printf("There is no option for flag %s\n",argv[i]);
	}
	else { // this is an arg
	    char* end;
            if(number_argument == 0) { //first number argument
		order = strtol(argv[i], &end, 10);
		if(*end != '\0' || order < 2) {
		    incorrect_usage_of_args = 1;
		    printf("Order is not a number bigger then 2!\n");
		}
	    }
	    else if(number_argument == 1) { //second number argument
		number_of_threads = strtol(argv[i], &end, 10);
		if(*end != '\0' || number_of_threads < 2) {
		    incorrect_usage_of_args = 1;
		    printf("Number of Threads is not a number bigger then 2!\n");
		}
	    }
	    else if(number_argument == 2) { //second number argument
		seed = strtol(argv[i], &end, 10);
		if(*end != '\0') {
		    incorrect_usage_of_args = 1;
		    printf("Seed is not a number!\n");
		}
	    }
	    else {
                incorrect_usage_of_args = 1;
		printf("Too many arguments!\n");
	    }
	    number_argument++;
	}
    }
    if(number_argument < 3) {
        printf("Arguments missing!\n");
	incorrect_usage_of_args = 1;
    }
    
    if(incorrect_usage_of_args == 1) {
        printf("Usage:\n");
        printf("$ ./jacobipar <N> <T> <seed>\n");
        printf("\tN - Matrix order (>= 2)\n");
        printf("\tT - Number of threads (>= 2)\n");
        printf("\tseed - seed for the pseudorandom number generator (>= 0)\n");
        printf("Flags:\n");
        printf("\t-h - displays help message\n");
        printf("\t-q - doesn't print matrix values nor result\n");
        exit(-1);
    }
    if(!quiet_mode)
        printf("Seed: %d\n", seed);

    srand(seed);
    
    matrix_t* A = init_rand_matrix(order, order);
    matrix_t* B = init_rand_matrix(order, 1);
    
    if(!quiet_mode) {
        print_matrix(A, 0);
        print_matrix(B, 0);
    }

    matrix_t* C = gaussjacobi_parallel(A, B, number_of_threads);
    
    if(!quiet_mode) {
        printf("Result: \n");
        print_matrix(C, 0);
    }

    free_matrix(C);
    free_matrix(A);
    free_matrix(B);

    return 0;
 }
