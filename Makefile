all:
	gcc src/jacobiseq.c src/matrix.c -I include/ -fopenmp -O3 -g -o jacobiseq
	gcc src/jacobipar.c src/matrix.c -I include/ -fopenmp -O3 -o jacobipar

clean:
	rm -rf jacobiseq
	rm -rf jacobipar
	rm -rf *.out
	rm -rf perf.*
	rm -rf gecko_profile.json 

perf:
	perf record -g -- ./jacobipar 10000 10 42
	perf script report gecko

run:
	./jacobiseq 3 42
