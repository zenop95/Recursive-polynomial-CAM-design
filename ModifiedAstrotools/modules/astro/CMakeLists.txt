add_shared_library(LIB_NAME astro
    LINK_LIBRARIES
        ${CSPICE_LIBRARY}
        ${DACE_LIBRARIES}
        geco
    INCLUDE_DIRS
        ${CSPICE_INCLUDE_DIR}
        ${DACE_INCLUDE_DIR}
    LIB_PUBLIC_HEADERS
        AstroCnvRef.h
        AstroDynModel.h
        AstroLibrary.h
        AstroRoutines.h
        Reference.h
    LIB_PRIVATE_SOURCES
        Reference.cpp
)
