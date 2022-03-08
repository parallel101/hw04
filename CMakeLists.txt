cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main main.cpp)

if (CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(main PRIVATE $<$<CONFIG:Release>:-march=native -funroll-loops -O3>)
endif()
if (MSVC)
    target_compile_options(main PRIVATE $<$<CONFIG:Release>:/arch:AVX2 /fp:fast>)
endif()
