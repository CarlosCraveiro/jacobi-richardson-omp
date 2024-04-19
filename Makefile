all:
	gcc src/jacobiseq.c src/matrix.c -I include/ -o jacobiseq

clean:
	rm -rf jacobiseq
	rm -rf jacobipar
