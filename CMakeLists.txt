cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

project(ParDeTria CXX)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()

# set output dir to bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

# setup ccache if found
find_program(CCACHE_PROG ccache)
if (NOT CCACHE_PROG MATCHES NOTFOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif ()

# setup modern c++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++17 -pedantic -Wall -Wextra -Werror")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c11 -pedantic -Wall -Wextra -Werror")

# add project include directories
include_directories("${PROJECT_SOURCE_DIR}/include")

# add all sources
add_subdirectory(src)
