set(cavity_SOURCES 
    GePolCavity.cpp 
    TsLessCavity.cpp
    WaveletCavity.cpp
    CavityFactory.cpp
)

add_library(cavity ${cavity_SOURCES})

install(TARGETS cavity ARCHIVE DESTINATION lib)
 
install(FILES Cavity.hpp DESTINATION include)
install(FILES GePolCavity.hpp DESTINATION include)
install(FILES TsLessCavity.hpp DESTINATION include)
install(FILES WaveletCavity.hpp DESTINATION include)
install(FILES CavityFactory.hpp DESTINATION include)
