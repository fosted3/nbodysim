#!/bin/sh
#This will pass your first bash arg to valgrind. Common usage is -v
make clean && make
valgrind --leak-check=full --show-reachable=yes --read-var-info=yes $1 ./bin/nbodysim
