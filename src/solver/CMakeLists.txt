set(solver_SOURCES
	CPCMSolver.cpp
  	IEFSolver.cpp
  	PWCSolver.cpp
  	PWLSolver.cpp
	SolverFactory.cpp
  	WEMSolver.cpp
)

if(ENABLE_TESTS)
	add_library(solver ${solver_SOURCES})
	add_dependencies(solver cavity green wem pwl)
	install(TARGETS solver ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${solver_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(FILES SolverFactory.hpp DESTINATION include)
install(FILES IEFSolver.hpp DESTINATION include)
install(FILES CPCMSolver.hpp DESTINATION include)
install(FILES PCMSolver.hpp DESTINATION include)
install(FILES PWCSolver.hpp DESTINATION include)
install(FILES PWLSolver.hpp DESTINATION include)
install(FILES WEMSolver.hpp DESTINATION include)
