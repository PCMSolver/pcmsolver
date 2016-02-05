# List of headers
list(APPEND headers_list CPCMSolver.hpp IEFSolver.hpp PCMSolver.hpp RegisterSolverToFactory.hpp)

# List of sources
list(APPEND sources_list CPCMSolver.cpp IEFSolver.cpp)

add_library(solver OBJECT ${sources_list} ${headers_list})
set_target_properties(solver PROPERTIES POSITION_INDEPENDENT_CODE 1 INTERPROCEDURAL_OPTIMIZATION 1)
set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/solver)
endforeach()
