cmake_minimum_required(VERSION 3.10)
project(acwv)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# find_package(MPI REQUIRED)

# set(CMAKE_CXX_COMPILER mpicxx)
set (CMAKE_CXX_COMPILER g++)

file(GLOB source_files src/main.cpp src/*2d.cpp)

include_directories(include)

add_executable(acwv ${source_files})
target_link_libraries(acwv ${VTK_LIBRARIES})
