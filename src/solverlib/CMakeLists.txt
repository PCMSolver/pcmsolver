set (solverlib_SOURCES
PCMSolver.cpp
)

add_library(solver ${solverlib_SOURCES})

install(TARGETS solver
	ARCHIVE DESTINATION lib)
