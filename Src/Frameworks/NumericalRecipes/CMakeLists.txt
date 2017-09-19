set(HEADERS nr3.h)
set(SOURCES nr3.cpp)

add_library(NR3 ${HEADERS} ${SOURCES})
target_include_directories(NR3 PUBLIC ./)
target_link_libraries(NR3 PUBLIC)
