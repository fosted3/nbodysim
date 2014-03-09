#!/bin/sh
g++ -g -o main *.cpp
valgrind --leak-check=full --show-reachable=yes --read-var-info=yes ./main
