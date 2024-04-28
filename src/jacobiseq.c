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

matrix_value_t gaussjacobi_error(const matrix_t* Xk, const matrix_t* Xkprev) {
    matrix_value_t *diff_vec = malloc(Xk->rows * sizeof(*diff_vec));
    matrix_value_t *abs_xk_vec = malloc(Xk->rows * sizeof(*abs_xk_vec));
    
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

matrix_t* gaussjacobi(const matrix_t* A, const matrix_t* B) {
    matrix_t* Xk = init_matrix(B->rows, 1, 1);
    matrix_t* Xkprev = init_matrix(B->rows, 1, 1);
    int itr = 0;
    do {
        matrix_swap(Xkprev, Xk); 
        
        for(int i = 0; i < B->rows; i++) {
            matrix_value_t xi = B->data[B->columns * i + 0];

            for(int j = 0; j < A->columns ; j++) {
                if(i != j) xi += -1 * A->data[A->columns * i + j] * Xkprev->data[Xkprev->columns * j + 0];
            }

            Xk->data[Xk->columns * i + 0] = xi / A->data[A->columns * i + i];
        }
        //printf("Iteration %d\n", itr++);
        //print_matrix(&Xk, 1);
    
    } while(gaussjacobi_error(Xk, Xkprev) > THRESHOLD);
    
    free_matrix(Xkprev);

    return Xk;
}

int main(int argc, char* argv[]) { 
    int order, seed;
    int quiet_mode = 0;
    int incorrect_usage_of_args = 0;
    int number_argument = 0;
    for(int i = 1; i < argc; i++) {
        if(argv[i][0] == '-') { // this is a flag
	    if(strcmp(argv[i],"-q") == 0) quiet_mode = 1; //quiet flag
	    else if(strcmp(argv[i],"-h") == 0) incorrect_usage_of_args = 1; //force help message
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
    if(number_argument < 2) {
        printf("Arguments missing!\n");
	incorrect_usage_of_args = 1;
    }
    
    if(incorrect_usage_of_args == 1) {
        printf("Correct Usage:\n");
        printf("$ ./jacobiseq <N> <seed>\n");
        printf("\tN - Matrix order (>= 2)\n");
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
    matrix_t* C = gaussjacobi(A, B);
    
    if(!quiet_mode) {
        printf("Result: \n");
        print_matrix(C, 0);
    }

    free_matrix(C);
    free_matrix(A);
    free_matrix(B);

    return 0;
 }
