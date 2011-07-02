set (solverlib_SOURCES PCMSolver.cpp)

add_executable(pcmsolver.x test.cpp)
target_link_libraries(pcmsolver.x solver metal green pedra cavity getkw)
target_link_libraries(pcmsolver.x ${LIBS})

find_package(PythonInterp REQUIRED)

set (INST_LIBDIR ${CMAKE_SOURCE_DIR}/external/share/libgetkw)
set (PCMSOLVER_EXECUTABLE ${CMAKE_BINARY_DIR}/pcmsolver.x)
Option (PCMSOLVER_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
#ask jonas about path
configure_file(${PROJECT_SOURCE_DIR}/src/solverlib/pcmsolver.in pcmsolver)
execute_process(COMMAND chmod 755 ${CMAKE_BINARY_DIR}/src/solverlib/pcmsolver OUTPUT_QUIET) 



add_library(solver ${solverlib_SOURCES})

install(TARGETS solver
	ARCHIVE DESTINATION lib)

install(TARGETS pcmsolver.x
	RUNTIME DESTINATION bin
	)

install(PROGRAMS pcmsolver 
	DESTINATION bin
	)

