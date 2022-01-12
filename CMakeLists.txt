cmake_minimum_required(VERSION 3.12)
project(hellocmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main main.cpp)

add_executable(opt_main opt_main.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(opt_main PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(opt_main PUBLIC -O3  -ffast-math -march=native)
# -march=native 让编译器自动检测当前硬件支持的指令集，比如avm的ymm