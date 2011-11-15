set (solverlib_SOURCES 
  IEFSolver.cpp 
  interface.cpp
  PCMSolver.cpp
  WEMSolver.cpp
)

add_executable(pcmsolver.x test.cpp)

target_link_libraries(pcmsolver.x solver green cavity pedra metal wem wavcav getkw)
target_link_libraries(pcmsolver.x ${LIBS})

find_package(PythonInterp REQUIRED)

set (INST_LIBDIR ${PROJECT_SOURCE_DIR}/external/share/libgetkw)
#set (PCMSOLVER_EXECUTABLE ${PROJECT_BINARY_DIR}/pcmsolver.x)
Option (PCMSOLVER_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
#ask jonas about path
configure_file(${PROJECT_SOURCE_DIR}/src/solverlib/pcmsolver.in pcmsolver)
execute_process(COMMAND chmod 755 src/solverlib/pcmsolver OUTPUT_QUIET) 

add_library(solver ${solverlib_SOURCES})

install(TARGETS solver
	ARCHIVE DESTINATION lib)

install(TARGETS pcmsolver.x
	RUNTIME DESTINATION bin
	)

install(PROGRAMS ${PROJECT_BINARY_DIR}/pcmsolver 
	DESTINATION bin
	)

