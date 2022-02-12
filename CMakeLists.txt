cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)

add_executable(main main.cpp)
target_compile_options(main PUBLIC -ffast-math -march=native)
target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
