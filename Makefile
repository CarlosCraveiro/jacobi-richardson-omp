#************************************************************
#          Trabalho Prático 1 - SCC0903                     *
#                                                           *
#      Nome: Artur Brenner Weber                            *
#      nUSP: 12675451                                       *
#      Nome: Carlos Henrique Craveiro Aquino Veras          *
#      nUSP: 12547187                                       *
#      Nome: Gabriel Franceschi Libardi                     *
#      nUSP: 11760739                                       *
#      Nome: Ivan Roberto Pancheniak                        *
#      nUSP: 12624224                                       *
#      Data de última atualizacao: 15/5/2024                *
#      Ambiente: VSCode 1.89.1                              *
#                                                           *
#           Arquivo Makefile                                *
#***********************************************************/

# Flags de Compilação

CFLAGS        = -Ofast -ftree-vectorize \
                -mtune=native -fopenmp-simd -fopt-info-optimized=stdout
LDFLAGS       = -fopenmp -lm
INCLUDE_PATHS = include/

all:
	gcc src/jacobiseq.c src/matrix.c $(CFLAGS) -I $(INCLUDE_PATHS) $(LDFLAGS) -o jacobiseq
	gcc src/jacobipar.c src/matrix.c $(CFLAGS) -I $(INCLUDE_PATHS) $(LDFLAGS) -o jacobipar

clean:
	rm -rf jacobiseq
	rm -rf jacobipar
	rm -rf *.out
	rm -rf perf.*
	rm -rf gecko_profile.json 

perf:
	perf record -g -- ./jacobipar 20000 10 42 5
	perf script report gecko

run:
	./jacobiseq 20000 42 0

parallel:
	./jacobipar 2000 10 42 5
