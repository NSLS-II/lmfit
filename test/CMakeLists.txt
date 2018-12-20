set(unit_tests
    test_qr
    test_inversion
    )

foreach(app ${unit_tests})
    add_executable(${app} ${app}.c)
    target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${lmfit_SOURCE_DIR}/lib)
    target_link_libraries(${app} ${lmfit_LIBRARY})
endforeach()
