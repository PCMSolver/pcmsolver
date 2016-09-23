# List of headers
list(APPEND headers_list CPCMSolver.hpp IEFSolver.hpp PCMSolver.hpp RegisterSolverToFactory.hpp SolverData.hpp)

# List of sources
list(APPEND sources_list CPCMSolver.cpp IEFSolver.cpp)

add_library(solver OBJECT ${sources_list} ${headers_list})
set_target_properties(solver PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(solver generate-config-hpp)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN}/solver)
endforeach()
