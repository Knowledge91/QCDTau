set(SOURCES CRunDec.cpp)
set(HEADERS CRunDec.h)

add_library(CRunDec ${SOURCES} ${HEADERS})
target_include_directories(CRunDec PUBLIC ./)
target_link_libraries(CRunDec PUBLIC)
