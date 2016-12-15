# List of headers
list(APPEND headers_list BIOperatorData.hpp BoundaryIntegralOperator.hpp IBoundaryIntegralOperator.hpp Collocation.hpp Purisima.hpp Numerical.hpp)

# List of sources
list(APPEND sources_list IBoundaryIntegralOperator.cpp Collocation.cpp Purisima.cpp Numerical.cpp)

add_library(bi_operators OBJECT ${headers_list} ${sources_list})
set_target_properties(bi_operators PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(bi_operators generate-config-hpp)
target_compile_options(bi_operators PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/bi_operators)
endforeach()
