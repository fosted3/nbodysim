#!/bin/sh
make clean && make
valgrind --leak-check=full --show-reachable=yes --read-var-info=yes -v ./bin/nbodysim
