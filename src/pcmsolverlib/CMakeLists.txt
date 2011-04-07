set (pcmsolverlib_SOURCES
PCMSolver.cpp
)

add_library(pcmsolver ${pcmsolverlib_SOURCES})

install(TARGETS pcmsolver
	ARCHIVE DESTINATION lib)
