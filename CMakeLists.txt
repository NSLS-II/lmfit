cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

project(lmfit C)
set(lmfit_SOVERSION                  8) # API version
set(lmfit_VERSION ${lmfit_SOVERSION}.1) # lib version

# --- Declare project-wide user flags, and set default values ---
option(FITTEST "Build with FitTest" OFF)
option(LIB_MAN "Build library man pages" ON)
option(PEDANTIC "Compile with pedantic warnings" ON)

set(destination ${CMAKE_INSTALL_PREFIX})

if(MSVC)
    add_compile_options(-O2)
    if(PEDANTIC)
        add_compile_options(-pedantic -Wall)
    endif()
endif()

configure_file("lmfit.pc.in" "lmfit.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lmfit.pc"
    DESTINATION "${destination}/lib/pkgconfig/")

add_subdirectory(lib)
add_subdirectory(demo)
if (${LIB_MAN})
    add_subdirectory(man)
endif()
if (${FITTEST})
    find_package(FitTest REQUIRED)
    add_subdirectory(fittest)
endif()
