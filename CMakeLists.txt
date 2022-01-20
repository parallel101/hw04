cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main main.cpp)
#5th: 设置-ffast-math有奇效
#13th: 设置-march=native 也有奇效. 至今为之，加速了10倍
target_compile_options(main PUBLIC -ffast-math -march=native -O3)

