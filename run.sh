!/bin/sh
set -e
cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=True
cmake --build build 
build/main
# build/opt_main

g++ -std=c++17 -march=native -ffast-math -O3 -fopenmp  opt_main.cpp -o opt_main
./opt_main