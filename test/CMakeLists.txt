# interactive programs for trying things during development
set(runnables
    run_qr
    )

foreach(app ${runnables})
    add_executable(${app} ${app}.c)
    target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${lmfit_SOURCE_DIR}/lib)
    target_link_libraries(${app} ${lmfit_LIBRARY})
endforeach()

# unit tests
set(unit_tests
    test_inversion
    )

foreach(app ${unit_tests})
    add_executable(${app} ${app}.c)
    add_test(${app} ${CMAKE_CURRENT_BINARY_DIR}/${app})
    target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${lmfit_SOURCE_DIR}/lib)
    target_link_libraries(${app} ${lmfit_LIBRARY})
endforeach()

# functional tests
add_test(curve1 ${CMAKE_BINARY_DIR}/demo/curve1)
add_test(surface1 ${CMAKE_BINARY_DIR}/demo/surface1)
