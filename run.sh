#!/bin/sh
set -e

g++ -std=c++17 -march=native -ffast-math -O3 -fopenmp  main.cpp -o opt_main
#g++ -std=c++17 -march=native -O3 -fopenmp  main.cpp -o opt_main
./opt_main

