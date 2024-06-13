cmake_minimum_required(VERSION 3.10)
project(acwv)

find_package(PkgConfig REQUIRED)

if(DEFINED ENV{PETSC_DIR})
  set (ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  pkg_search_module(PETSc REQUIRED PETSc>3.10)
  if(PETSc_FOUND)
     message(STATUS "PETSc ${PETSc_VERSION} was found!")
     message(STATUS "PETSc_INCLUDE: ${PETSc_INCLUDE_DIRS} ")
     message(STATUS "PETSc_LIBRARIES: ${PETSc_LINK_LIBRARIES} ")
     SET(PETSC_ENABLE 1)
     include_directories(${PETSc_INCLUDE_DIRS})
     link_directories(${PETSc_LIBRARY_DIRS})
  endif()
else()
  message(STATUS "PETSc was not founded! Please define PETSC_DIR var!!")
endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

file(GLOB source_files src/main.cpp src/*2d.cpp)

include_directories(include)

add_executable(acwv ${source_files})
target_link_libraries(acwv ${VTK_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSc_LIBRARIES})