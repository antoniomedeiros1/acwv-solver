cmake_minimum_required(VERSION 3.10)
project(acwv)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

find_package(PETSc REQUIRED)
include_directories(${PETSc_INCLUDE_DIRS})
link_directories(${PETSc_LIBRARY_DIRS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

file(GLOB source_files src/main.cpp src/*2d.cpp)

include_directories(include)

add_executable(acwv ${source_files})
target_link_libraries(acwv ${VTK_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSc_LIBRARIES})