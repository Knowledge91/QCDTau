################################
# Testing
################################
# has 'project(gtest)'.
add_subdirectory(googletest)

enable_testing()

# Include the gtest library. gtest_SOURCE_DIR is available due to
# 'project(gtest)' above.
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

##############
# Unit Tests
##############
file(GLOB SOURCES ../src/*.cpp ./*.cpp)

add_executable(runUnitTests ${SOURCES})
target_link_libraries(runUnitTests gtest gtest_main NR3 ALEPH CRunDec ${Boost_LIBRARIES} ${ROOT_LIBRARIES})
