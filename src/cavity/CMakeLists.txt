set(cavity_SOURCES 
    Cavity.cpp
    GePolCavity.cpp
    WaveletCavity.cpp
    CavityFactory.cpp
)

add_library(cavity ${cavity_SOURCES})
#add_dependencies(cavity external)

install(TARGETS cavity ARCHIVE DESTINATION lib)
 
install(FILES Cavity.hpp DESTINATION include)
install(FILES GePolCavity.hpp DESTINATION include)
install(FILES WaveletCavity.hpp DESTINATION include)
install(FILES CavityFactory.hpp DESTINATION include)
