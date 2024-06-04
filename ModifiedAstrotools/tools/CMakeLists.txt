add_tool(TOOL_NAME polyProp
    LINK_LIBRARIES
        ${DACE_LIBRARIES}
        ${JSON_LIBRARIES}
        atmos
        astro
        cfg
        dynorb
        geco
        odeprop
    TOOL_SOURCES
        polyProp.cpp
    )

add_tool(TOOL_NAME CGTProp
    LINK_LIBRARIES
        ${DACE_LIBRARIES}
        ${JSON_LIBRARIES}
        atmos
        astro
        cfg
        dynorb
        geco
        odeprop
    TOOL_SOURCES
        CGTProp.cpp
    )

add_tool(TOOL_NAME validatePoly
    LINK_LIBRARIES
        ${DACE_LIBRARIES}
        ${JSON_LIBRARIES}
        atmos
        astro
        cfg
        dynorb
        geco
        odeprop
    TOOL_SOURCES
        validatePoly.cpp
    )