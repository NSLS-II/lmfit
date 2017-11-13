cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(lmfit)

set(lmfit_VERSION 3.5)

set(destination_lib ${CMAKE_INSTALL_PREFIX}/lib)

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -pedantic -Wall -Wno-sign-compare -Wno-unused-result -Wno-parentheses  -Wno-unknown-pragmas -Werror")

add_subdirectory(lib)
add_subdirectory(demo)
