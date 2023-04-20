cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main main.cpp)

if(MSVC)
else()
target_compile_options(main PUBLIC
    -Wall
    -Werror
    -O3
    -march=native
    -ffast-math
)
endif()