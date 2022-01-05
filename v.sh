set -e

g++ -march=native -O3 -fopenmp -fomit-frame-pointer -fverbose-asm -S opt_main.cpp -o /tmp/opt_main.S
# vim /tmp/opt_main.S
cat /tmp/opt_main.S > a.txt