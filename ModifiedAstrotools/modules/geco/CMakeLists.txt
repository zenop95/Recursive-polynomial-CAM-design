add_shared_library(LIB_NAME geco
    LINK_LIBRARIES
        ${CSPICE_LIBRARY}
        ${JSON_LIBRARY}
    INCLUDE_DIRS
        ${JSON_INCLUDE_DIR}
    LIB_PUBLIC_HEADERS
        Factory.hpp
        JsonParser.hpp
        LazyLoader.hpp
        MapEnum.h
        NonCopyable.hpp
    LIB_PRIVATE_HEADERS
        LazyLoader_impl.hpp
    LIB_SOURCES
        JsonParser.cpp
        MapEnum.cpp
)