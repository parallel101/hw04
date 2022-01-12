cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(main PRIVATE -march=native -ffast-math -fverbose-asm -fomit-frame-pointer)
# target_compile_options(main PRIVATE -march=native -fverbose-asm -fomit-frame-pointer)