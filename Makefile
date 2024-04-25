all:
	gcc src/jacobiseq.c src/matrix.c -I include/ -fopenmp -O3 -o jacobiseq
	gcc src/jacobipar.c src/matrix.c -I include/ -fopenmp -O3 -o jacobipar

clean:
	rm -rf jacobiseq
	rm -rf jacobipar
	rm -rf *.out

run:
	./jacobiseq 3 42
