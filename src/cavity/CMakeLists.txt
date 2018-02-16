# List of headers
list(APPEND headers_list
Cavity.hpp
CavityData.hpp
Element.hpp
GePolCavity.hpp
ICavity.hpp
RestartCavity.hpp
)

# List of sources
list(APPEND sources_list
Element.cpp
GePolCavity.cpp
ICavity.cpp
RestartCavity.cpp
)

add_library(cavity OBJECT ${sources_list} ${headers_list})
set_target_properties(cavity
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
target_compile_definitions(cavity
  PUBLIC
    PCMSolver_EXPORTS
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(cavity custom_boost)
endif()
add_dependencies(cavity generate-config-hpp)
target_compile_options(cavity
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/cavity)
endforeach()
