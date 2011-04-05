set (pcmsolverlib_SOURCES
TraditionalPCMOperator.cpp)

add_library(pcmsolver ${pcmsolverlib_SOURCES})

install(TARGETS pcmsolver
	ARCHIVE DESTINATION lib)
