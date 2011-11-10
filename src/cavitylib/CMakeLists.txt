set(cavity_SOURCES 
    Cavity.cpp
    GePolCavity.cpp
	WaveletCavity.cpp
)

add_library(cavity ${cavity_SOURCES})

install(TARGETS cavity
	ARCHIVE DESTINATION lib)

install(FILES Cavity.h DESTINATION include)
install(FILES GePolCavity.h DESTINATION include)
install(FILES WaveletCavity.h DESTINATION include)
