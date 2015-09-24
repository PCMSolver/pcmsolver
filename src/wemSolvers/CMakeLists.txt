# List of headers
list(APPEND headers_list Cubature.hpp GaussLegendre.hpp GaussSquare.hpp Quadrature.hpp WEMPCG.hpp WEMPGMRES.hpp WEMRHS.hpp)

# List of sources
list(APPEND sources_list GaussLegendre.cpp GaussSquare.cpp WEMPCG.cpp WEMPGMRES.cpp WEMRHS.cpp)

# Update bar chart
if(BUILD_DOCUMENTATION)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR})
endif()

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/wemSolvers)
endforeach()
