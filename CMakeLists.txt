cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main main.cpp)
# target_compile_options(main PUBLIC -ffast-math -march=native)
add_executable(work work.cpp)
target_compile_options(work PUBLIC -ffast-math -march=native)

add_executable(work2 work2.cpp)
target_compile_options(work2 PUBLIC -ffast-math -march=native)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
