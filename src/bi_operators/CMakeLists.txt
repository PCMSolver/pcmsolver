# List of headers
list(APPEND headers_list
BIOperatorData.hpp
BoundaryIntegralOperator.hpp
Collocation.hpp
IBoundaryIntegralOperator.hpp
Numerical.hpp
Purisima.hpp
)

# List of sources
list(APPEND sources_list
Collocation.cpp
IBoundaryIntegralOperator.cpp
Numerical.cpp
Purisima.cpp
)

add_library(bi_operators OBJECT ${sources_list} ${headers_list})
set_target_properties(bi_operators
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
add_dependencies(bi_operators generate-config-hpp)
target_compile_options(bi_operators
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/bi_operators)
endforeach()
