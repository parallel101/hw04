cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG"-g")
set(CMAKE_CXX_FLAGS_RELEASE"-O3")

# target_compile_options(testbench PUBLIC -ffast-math -march=native)
# find_package(OpenMP REQUIRED)
# target_link_libraries(testbench PUBLIC OpenMP::OpenMP_CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

    add_executable(main main.cpp)
