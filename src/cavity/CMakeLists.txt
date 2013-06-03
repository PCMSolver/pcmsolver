set(cavity_SOURCES 
    Cavity.cpp
    GePolCavity.cpp
    WaveletCavity.cpp
    SurfaceFunction.cpp
    CavityFactory.cpp
)

add_library(cavity ${cavity_SOURCES})
add_dependencies(cavity external)

install(TARGETS cavity
 	ARCHIVE DESTINATION lib)
 
install(FILES Cavity.h DESTINATION include)
install(FILES GePolCavity.h DESTINATION include)
install(FILES WaveletCavity.h DESTINATION include)
install(FILES SurfaceFunction.h DESTINATION include)
install(FILES CavityFactory.h DESTINATION include)
