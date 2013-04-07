set(solver_SOURCES 
  IEFSolver.cpp
  CPCMSolver.cpp 
  Interface.cpp
  PCMSolver.cpp
  PWCSolver.cpp
  PWLSolver.cpp
  WEMSolver.cpp
)

#add_executable(pcmsolver.x test.cpp)
#
#target_link_libraries(pcmsolver.x solver green cavity pedra metal wem wavcav getkw)
#target_link_libraries(pcmsolver.x ${LIBS})



add_library(solver ${solver_SOURCES})
add_dependencies(solver external)

install(TARGETS solver
	ARCHIVE DESTINATION lib)

install(FILES IEFSolver.h DESTINATION include)
install(FILES CPCMSolver.h DESTINATION include)
install(FILES Interface.h DESTINATION include)
install(FILES PCMSolver.h DESTINATION include)
install(FILES PWCSolver.h DESTINATION include)
install(FILES PWLSolver.h DESTINATION include)
install(FILES WEMSolver.h DESTINATION include)
