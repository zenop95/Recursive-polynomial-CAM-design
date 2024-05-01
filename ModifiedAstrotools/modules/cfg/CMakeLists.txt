add_shared_library(LIB_NAME cfg
    LINK_LIBRARIES
        ${CSPICE_LIBRARY}
        ${JSON_LIBRARY}
        geco
    INCLUDE_DIRS
        ${JSON_INCLUDE_DIR}
    LIB_PUBLIC_HEADERS
        Config.hpp
        Environment.hpp
        EnvironmentOption.hpp
    LIB_PRIVATE_HEADERS
        EnvironmentOption_impl.hpp
    LIB_SOURCES
        Config.cpp
        Environment.cpp
)