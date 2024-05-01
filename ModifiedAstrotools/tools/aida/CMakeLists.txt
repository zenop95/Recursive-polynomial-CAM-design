add_tool(TOOL_NAME runAIDA
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
        RunAida.cpp
    )