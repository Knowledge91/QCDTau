set(SOURCES aleph14_vpa.f90 num_const.f90)

add_library(ALEPH ${SOURCES})
target_link_libraries(ALEPH PUBLIC)
