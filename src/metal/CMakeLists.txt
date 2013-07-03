set(metal_SOURCES
    gsfera.F
    plgndr.F
   )

add_library(metal ${metal_SOURCES})

install(TARGETS metal ARCHIVE DESTINATION lib)

