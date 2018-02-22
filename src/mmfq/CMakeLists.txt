# List of headers
list(APPEND headers_list
  FQOhno.hpp
  )

# List of sources
list(APPEND sources_list
  FQOhno.cpp
  )

add_library(mmfq OBJECT ${sources_list} ${headers_list})
set_target_properties(mmfq
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
target_compile_definitions(mmfq
  PUBLIC
    PCMSolver_EXPORTS
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(mmfq custom_boost)
endif()
add_dependencies(mmfq generate-config-hpp)
target_compile_options(mmfq
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/mmfq)
endforeach()
