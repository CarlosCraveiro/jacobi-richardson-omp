#!/bin/sh
# from https://stackoverflow.com/a/69430652 
for i in {1..30}; do time $@; done 2>&1 | grep ^real | sed s/,/./ | sed -e s/.*m// | awk '{sum += $1} END {print sum / NR}'
