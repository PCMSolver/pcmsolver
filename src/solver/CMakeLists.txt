# List of headers
list(APPEND headers_list CPCMSolver.hpp IEFSolver.hpp PCMSolver.hpp RegisterSolverToFactory.hpp)

# List of sources
list(APPEND sources_list CPCMSolver.cpp IEFSolver.cpp)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/solver)
endforeach()
