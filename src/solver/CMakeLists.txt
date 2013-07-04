set(solver_SOURCES
	SolverFactory.cpp
  	IEFSolver.cpp
	CPCMSolver.cpp 
  	PCMSolver.cpp
  	PWCSolver.cpp
  	PWLSolver.cpp
  	WEMSolver.cpp
)


add_library(solver ${solver_SOURCES})

install(TARGETS solver ARCHIVE DESTINATION lib)

install(FILES SolverFactory.h DESTINATION include)
install(FILES IEFSolver.h DESTINATION include)
install(FILES CPCMSolver.h DESTINATION include)
install(FILES PCMSolver.h DESTINATION include)
install(FILES PWCSolver.h DESTINATION include)
install(FILES PWLSolver.h DESTINATION include)
install(FILES WEMSolver.h DESTINATION include)
