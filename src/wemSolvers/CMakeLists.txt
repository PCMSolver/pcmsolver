# List of headers
list(APPEND headers_list Cubature.hpp GaussLegendre.hpp GaussSquare.hpp Quadrature.hpp WEMPCG.hpp WEMPGMRES.hpp WEMRHS.hpp)

# List of sources
list(APPEND sources_list GaussLegendre.cpp GaussSquare.cpp WEMPCG.cpp WEMPGMRES.cpp WEMRHS.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/wemSolvers)
endforeach()
