cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -fno-optimize-sibling-calls -ffast-math -march=native -fopt-info-vec-all -fopt-info-vec-missed -fdump-tree-vect-all  -mssse3 -msse4.1 -msse4.2")
# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -fno-optimize-sibling-calls -ffast-math -march=native -fopt-info-vec-all -fopt-info-vec-missed -fdump-tree-vect-all  -mssse3 -msse4.1 -msse4.2")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math  -mavx512f")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math  -mavx512f")

add_executable(main main.cpp)
