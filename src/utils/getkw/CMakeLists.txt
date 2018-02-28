# List of headers
list(APPEND _headers_list
Getkw.h
Section.h
Keyword.h
GetkwError.h
)

# List of sources
list(APPEND _sources_list
Getkw.cpp
GetkwError.cpp
Section.cpp
messages.cpp
)

add_library(getkw OBJECT ${_sources_list} ${_headers_list})
set_target_properties(getkw
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
target_compile_definitions(getkw
  PUBLIC
    PCMSolver_EXPORTS
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(getkw custom_boost)
endif()
add_dependencies(getkw generate-config-hpp)
target_compile_options(getkw
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )

# Sets install directory for all the headers in the list
foreach(_header ${_headers_list})
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils/getkw)
endforeach()
