include_directories(${CMAKE_CURRENT_BINARY_DIR} ${lmfit_SOURCE_DIR}/lib)

set(demos
    curve1
    fourexp
    hat
    morobropro
    powell
    surface1
    )

foreach(demo ${demos})
    add_executable(${demo} ${demo}.c)
    target_link_libraries(${demo} ${lmfit_LIBRARY})
endforeach()
