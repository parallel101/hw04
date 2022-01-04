cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 CONFIG REQUIRED)

add_executable(main main.cpp)
target_compile_options(main PUBLIC /arch:AVX2 /fp:fast /favor:INTEL64 /openmp)
target_link_libraries(main PRIVATE Eigen3::Eigen)
