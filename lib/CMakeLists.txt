include(GNUInstallDirs)
set(lmfit_LIBRARY lmfit PARENT_SCOPE)

set(HEADER_FILES lmcurve.h lmmin.h lmstruct.h)
set(SOURCE_FILES lmcurve.c lmmin.c lminvert.c)

add_library(lmfit ${SOURCE_FILES} ${HEADER_FILES})

set_target_properties(
    lmfit PROPERTIES
    OUTPUT_NAME lmfit
	VERSION ${lmfit_VERSION}
	SOVERSION ${lmfit_VERSION_MAJOR}
	PUBLIC_HEADER "${HEADER_FILES}"
)

target_include_directories(lmfit PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(LinkLibMath)
link_libm(lmfit)

if (LIB_INSTALL)
	include(CMakePackageConfigHelpers)

    install(
		TARGETS lmfit 
		EXPORT lmfit-targets
		LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
		PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT Libraries)

	install(
		EXPORT lmfit-targets
		FILE lmfit-targets.cmake
		NAMESPACE lmfit::
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lmfit"
	)

	configure_package_config_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/lmfit-config.cmake.in"
		"${CMAKE_CURRENT_BINARY_DIR}/lmfit-config.cmake"
		INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lmfit"
	)
	
	write_basic_package_version_file(
		"${CMAKE_CURRENT_BINARY_DIR}/lmfit-config-version.cmake"
		VERSION ${lmfit_VERSION}
		COMPATIBILITY AnyNewerVersion
	)

	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/lmfit-config.cmake"
		"${CMAKE_CURRENT_BINARY_DIR}/lmfit-config-version.cmake"
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lmfit"
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
