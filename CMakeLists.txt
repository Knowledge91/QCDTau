cmake_minimum_required (VERSION 3.9.2)
project (QCDTau CXX Fortran)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")

include_directories(SRC)
add_subdirectory(SRC)

add_executable(QCDTau main.cpp)
target_link_libraries(QCDTau PUBLIC CRunDec ALEPH)
