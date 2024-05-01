add_tool(TOOL_NAME dals
    LINK_LIBRARIES
        ${DACE_LIBRARIES}
        ${JSON_LIBRARIES}
        atmos
        astro
        dynorb
        odeprop
        dlib::dlib
        Eigen3::Eigen
    INCLUDE_DIRS
        ${DACE_INCLUDE_DIR}
        ${JSON_INCLUDE_DIR}
        ${CSPICE_INCLUDE_DIR}
    TOOL_SOURCES
        Main_DALS.cpp
    )