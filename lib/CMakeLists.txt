set(lib lmfit)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)
set(${lib}_LIBRARY_TYPE SHARED)

set(src_files lmcurve.c lmmin.c lminvert.c)
set(inc_files lmcurve.h lmmin.h lmstruct.h)

add_library(${lib} SHARED ${src_files})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${lmfit_VERSION_MAJOR}
    SOVERSION ${lmfit_VERSION_MINOR})

include(LinkLibMath)
link_libm(${lib})

if (LIB_INSTALL)
    install(
        TARGETS ${lib} LIBRARY
        DESTINATION ${destination}/lib
        COMPONENT Libraries)
    install(
        FILES ${inc_files}
        DESTINATION ${destination}/include
        COMPONENT Headers)
endif()
