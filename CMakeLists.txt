cmake_minimum_required (VERSION 3.9.2)

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." OFF) # Makes boolean 'test' available

project (QCDTau CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")

find_package(Boost 1.65.1)
find_package(ROOT REQUIRED COMPONENTS MINUIT)
# message(STATUS "${ROOT_USE_FILE}" )

# includes somehow the Root header files
include(${ROOT_USE_FILE})


################################
# Normal Libraries & Executables
################################
include_directories(SRC src/Frameworks ${Boost_INCLUDE_DIR})
add_subdirectory(SRC)

file(GLOB SOURCES src/*.cpp)

add_executable(QCDTau main.cpp ${SOURCES})
target_link_libraries(QCDTau PUBLIC CRunDec ALEPH NR3 ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

################################
# Testing
################################
if(test)
    add_subdirectory(test)
endif()
