# Jacobi-Richardson em OpenMP
Implementação do método de Jacobi-Richardson, um método iterativo para resolver sistema linear de equações, em OpenMP. Para compilar, basta executar o comando:
```
make all
```
Para executar o código sequencial, basta executar:
```
./jacobiseq <N> <seed> <row num>
```
Em que:
- \<N\> é a ordem da matriz.
- \<seed\> é a semente passada ao srand() para a geração de números pseudoaleatórios.
- \<row num\> é a linha na qual o termo constante calculado pelo método de gauss-jacobi é comparado com o seu valor real.

Para executar o código paralelo, basta rodar na linha de comando:
```
./jacobipar <N> <T> <seed> <row num>
```
Em que:
- \<N\> é a ordem da matriz.
- \<T\> é o número de threads.
- \<seed\> é a semente passada ao srand() para a geração de números pseudoaleatórios.
- \<row num\> é a linha na qual o termo constante calculado pelo método de gauss-jacobi é comparado com o seu valor real.