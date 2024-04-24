#!/bin/bash

./jacobiseq 3 42 > seq.out
./jacobipar 3 12 42 > par.out
diff seq.out par.out
