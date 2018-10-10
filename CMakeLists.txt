cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(lmfit C)
set(lmfit_SOVERSION                  8) # API version
set(lmfit_VERSION ${lmfit_SOVERSION}.0) # lib version

# --- Declare project-wide user flags, and set default values ---
option(FITTEST "Build with FitTest" OFF)
option(LIB_MAN "Build library man pages" ON)

set(destination ${CMAKE_INSTALL_PREFIX})

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -pedantic -Wall -Wno-sign-compare -Wno-unused-result -Wno-parentheses  -Wno-unknown-pragmas")

add_subdirectory(lib)
add_subdirectory(demo)
if (${LIB_MAN})
    add_subdirectory(man)
endif()
if (${FITTEST})
    find_package(FitTest REQUIRED)
    add_subdirectory(fittest)
endif()
