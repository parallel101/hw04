cmake_minimum_required(VERSION 3.12)
project(benchtest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


add_executable(main main.cpp)
find_package(OpenMP REQUIRED)
target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(main PUBLIC -ffast-math -march=native)

add_executable(baseline baseline.cpp)
target_link_libraries(baseline PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(baseline PUBLIC -ffast-math -march=native)
