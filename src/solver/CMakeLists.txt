set(solver_SOURCES
	SolverFactory.cpp
  	IEFSolver.cpp
	CPCMSolver.cpp 
  	PWCSolver.cpp
  	PWLSolver.cpp
  	WEMSolver.cpp
)


add_library(solver ${solver_SOURCES})

install(TARGETS solver ARCHIVE DESTINATION lib)

install(FILES SolverFactory.hpp DESTINATION include)
install(FILES IEFSolver.hpp DESTINATION include)
install(FILES CPCMSolver.hpp DESTINATION include)
install(FILES PCMSolver.hpp DESTINATION include)
install(FILES PWCSolver.hpp DESTINATION include)
install(FILES PWLSolver.hpp DESTINATION include)
install(FILES WEMSolver.hpp DESTINATION include)
