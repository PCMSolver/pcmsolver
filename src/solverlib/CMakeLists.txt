set(solver_SOURCES 
  IEFSolver.cpp 
  Interface.cpp
  PCMSolver.cpp
  PWCSolver.cpp
  PWLSolver.cpp
  Solvent.cpp
  WEMSolver.cpp
)

#add_executable(pcmsolver.x test.cpp)
#
#target_link_libraries(pcmsolver.x solver green cavity pedra metal wem wavcav getkw)
#target_link_libraries(pcmsolver.x ${LIBS})


set (INST_LIBDIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/share/libgetkw)
#set (PCMSOLVER_EXECUTABLE pcmsolver.x)
#set (PCMSOLVER_ERROR_LOG ${PROJECT_BINARY_DIR}/pcmsolver.err)
set (PCMSOLVER_EXECUTABLE ${PROJECT_BINARY_DIR}/pcmsolver.x)
Option (PCMSOLVER_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
#ask jonas about path
configure_file(pcmsolver.in ${PROJECT_BINARY_DIR}/bin/pcmsolver)
execute_process(COMMAND chmod 755 ${PROJECT_BINARY_DIR}/bin/pcmsolver OUTPUT_QUIET) 

add_library(solver ${solver_SOURCES})

install(TARGETS solver
	ARCHIVE DESTINATION lib)

install(FILES IEFSolver.h DESTINATION include)
install(FILES Interface.h DESTINATION include)
install(FILES PCMSolver.h DESTINATION include)
install(FILES PWCSolver.h DESTINATION include)
install(FILES PWLSolver.h DESTINATION include)
install(FILES Solvent.h   DESTINATION include)
install(FILES WEMSolver.h DESTINATION include)
