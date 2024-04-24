#!/bin/bash

# Loop de 0 a 100
for i in {0..100}
do
    for j in {0..100}
    do
        ./jacobiseq $i $j
    done
done
