# List of sources
target_sources(pcm-objlib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Collocation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IBoundaryIntegralOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Numerical.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Purisima.cpp
  )

# List of headers
list(APPEND headers_list
  BIOperatorData.hpp
  BoundaryIntegralOperator.hpp
  Collocation.hpp
  IBoundaryIntegralOperator.hpp
  Numerical.hpp
  Purisima.hpp
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS headers_list)
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/bi_operators)
endforeach()
