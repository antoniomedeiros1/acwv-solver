cmake_minimum_required(VERSION 3.10)
project(acwv)

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)

set(GLOB source_files
    src/main.cpp
    src/*2d.cpp
)

include_directories(include)

add_executable(acwv ${source_files})
