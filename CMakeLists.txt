cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

project(ParDeTria CXX)

# prohibit in-source builds
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif ()

# set output dir to bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

# setup ccache if found
find_program(CCACHE_PROG ccache)
if (NOT CCACHE_PROG MATCHES NOTFOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif ()

# adds the project-specific cmake module directory cmake/Modules to the cmake
# search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

# set default for MARCH
if(NOT MARCH)
    set(MARCH native)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSTIONS ON)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSTIONS ON)

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

# setup modern c++ flags
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${MARCH} -pedantic -Wall -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=${MARCH} -pedantic -Wall -Wextra")

# add project include directories
include_directories("${PROJECT_SOURCE_DIR}/includes")

option(VTUNE "Enable VTune instrumentation if found" OFF)
if (VTUNE)
    find_package(VTune)
    if (VTUNE_FOUND)
        message(STATUS "VTUNE instrumentation enabled")
        include_directories(
                SYSTEM
                ${VTUNE_INCLUDE_DIRS}
        )
        add_definitions(-DHAS_VTUNE)
    endif ()
endif ()

option(ADVISOR "Enable Advisor instrumentation if found" OFF)
if (ADVISOR)
    find_package(Advisor)
    if (Advisor_FOUND)
        message(STATUS "Advisor instrumentation enabled")
        include_directories(
                SYSTEM
                ${Advisor_INCLUDE_DIRS}
        )
        add_definitions(-DHAS_ADVISOR)
    endif ()
endif ()

find_package(Vc)
if (Vc_FOUND)
    message(STATUS "VC Library found")
    include_directories(
            SYSTEM
            ${Vc_INCLUDE_DIR}
    )
    add_definitions(-DHAS_Vc)
endif ()

option(TESTS "Enable tests" OFF)
if (TESTS)
    find_package(GTest)
    if(GTest_FOUND)
        message(STATUS "GTest enabled")
        include_directories(
                SYSTEM
                ${GTEST_INCLUDE_DIRS}
        )
    endif()
endif ()



# add all sources
add_subdirectory(src)
add_subdirectory(tests)

unset(MARCH CACHE)
