CFLAGS        = -Ofast -mprefer-vector-width=512 -ftree-vectorize -march=native\
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
