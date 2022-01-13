echo ===========Before===========
clang++ -std=c++17 -Wall -Wextra -Wpedantic -O3 main_before.cpp -o main_before
./main_before
echo ============After===========
clang++ -std=c++17 -Wall -Wextra -Wpedantic -O3 -ffast-math main.cpp -o main
./main

clang++ -std=c++17 -fomit-frame-pointer -fverbose-asm -S -O3 -ffast-math main.cpp -o main.s
