add_interface_library(LIB_NAME sada
    LINK_LIBRARIES
        ${CSPICE_LIBRARY}
        ${DACE_LIBRARIES}
        astro
    INCLUDE_DIRS
        ${CSPICE_INCLUDE_DIR}
        ${DACE_INCLUDE_DIR}
    LIB_PUBLIC_HEADERS
        constants.h
		Conversions.h
		HEOSATDynamics.h
		luniSolarShortPeriodics.h
		mean2oscNonSingular.h
		Transformations.h
	LIB_PRIVATE_HEADERS
        HEOSATDynamics_impl.h
)
