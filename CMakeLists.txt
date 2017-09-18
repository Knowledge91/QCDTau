cmake_minimum_required (VERSION 3.9.2)

project (QCDTau)
set(CMAKE_CXX_STANDARD 11)

include_directories(SRC)
add_subdirectory(SRC)

add_executable(QCDTau main.cpp)
target_link_libraries(QCDTau PUBLIC CRunDec)
