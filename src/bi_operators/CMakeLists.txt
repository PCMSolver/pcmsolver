# List of headers
list(APPEND headers_list BIOperatorData.hpp BoundaryIntegralOperator.hpp Collocation.hpp Purisima.hpp Numerical.hpp)

# List of sources
list(APPEND sources_list BoundaryIntegralOperator.cpp Collocation.cpp Purisima.cpp Numerical.cpp)

add_library(bi_operators OBJECT ${headers_list} ${sources_list})
set_target_properties(bi_operators PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(bi_operators generate-config-hpp)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/bi_operators)
endforeach()
