cmake_minimum_required (VERSION 3.9.2)

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." OFF) # Makes boolean 'test' available

project (QCDTau CXX Fortran)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")

find_package(Boost 1.65.1)
message(STATUS "${Boost_INCLUDE_DIR}" )

################################
# Normal Libraries & Executables
################################
include_directories(SRC Src/Frameworks ${Boost_INCLUDE_DIR} )
add_subdirectory(SRC)


add_executable(QCDTau main.cpp)
target_link_libraries(QCDTau PUBLIC CRunDec ALEPH NR3 ${Boost_LIBRARIES} )

################################
# Testing
################################
if(test)
    add_subdirectory(Test)
endif()
