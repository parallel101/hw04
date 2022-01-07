cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main main.cpp)
# target_compile_options(main PUBLIC -fopenmp -ffast-math -march=native)
add_executable(work work.cpp)
target_compile_options(work PUBLIC -fopenmp -ffast-math -march=native)
# target_compile_options(work PUBLIC -fopenmp -march=native)

find_package(Eigen3 REQUIRED)
add_executable(work_eigen work_eigen.cpp)
target_include_directories(work_eigen PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_options(work_eigen PUBLIC -march=native)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
