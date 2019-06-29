set(lib lmfit)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

set(src_files lmcurve.c lmmin.c lminvert.c)
set(inc_files lmcurve.h lmmin.h lmstruct.h)

add_library(${lib} ${src_files})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${lmfit_VERSION_MAJOR}
    SOVERSION ${lmfit_VERSION_MINOR})

include(LinkLibMath)
link_libm(${lib})

if (LIB_INSTALL)
	include(GNUInstallDirs)
	include(CMakePackageConfigHelpers)

    install(TARGETS ${lib} 
		EXPORT lmfit-targets
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT Libraries)

    install(FILES ${inc_files}
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Headers)

	install(EXPORT lmfit-targets
		FILE lmfit-targets.cmake
		NAMESPACE lmfit
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lmfit
	)

	configure_package_config_file(
		${CMAKE_CURRENT_SOURCE_DIR}/lmfit-config.cmake.in
		${CMAKE_CURRENT_BINARY_DIR}/lmfit-config.cmake
		INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lmfit
	)
	
	write_basic_package_version_file(
		${CMAKE_CURRENT_BINARY_DIR}/lmfit-config-version.cmake
		VERSION ${lmfit_VERSION_MAJOR}
		COMPATIBILITY AnyNewerVersion
	)

	install(
		FILES ${CMAKE_CURRENT_BINARY_DIR}/lmfit-config.cmake
		${CMAKE_CURRENT_BINARY_DIR}/lmfit-config-version.cmake
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lmfit
	)

    configure_file(
		"${CMAKE_SOURCE_DIR}/lmfit.pc.in" 
		"${CMAKE_CURRENT_BINARY_DIR}/lmfit.pc" @ONLY
	)


    install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/lmfit.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/"
	)

endif()
