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


set (INST_LIBDIR ${PROJECT_BINARY_DIR}/external/share/libgetkw)
#set (PCMSOLVER_EXECUTABLE pcmsolver.x)
#set (PCMSOLVER_ERROR_LOG ${PROJECT_BINARY_DIR}/pcmsolver.err)
set (PCMSOLVER_EXECUTABLE ${PROJECT_BINARY_DIR}/pcmsolver.x)
Option (PCMSOLVER_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
#ask jonas about path
configure_file(${PROJECT_SOURCE_DIR}/src/solverlib/pcmsolver.in pcmsolver)

add_library(solver ${solver_SOURCES})

install(TARGETS solver
	ARCHIVE DESTINATION lib)

install(FILES pcmsolver DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)

install(FILES IEFSolver.h DESTINATION include)
install(FILES Interface.h DESTINATION include)
install(FILES PCMSolver.h DESTINATION include)
install(FILES PWCSolver.h DESTINATION include)
install(FILES PWLSolver.h DESTINATION include)
install(FILES Solvent.h DESTINATION include)
install(FILES WEMSolver.h DESTINATION include)
