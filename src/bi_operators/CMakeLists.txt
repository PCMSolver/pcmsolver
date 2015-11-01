# List of headers
list(APPEND headers_list CollocationIntegrator.hpp IntegratorForward.hpp IntegratorHelperFunctions.hpp IntegratorTypes.hpp NumericalIntegrator.hpp PurisimaIntegrator.hpp)

# List of sources
list(APPEND sources_list )

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/bi_operators)
endforeach()
