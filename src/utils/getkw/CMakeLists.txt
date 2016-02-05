# List of headers
list(APPEND _headers_list Getkw.h Section.h Keyword.h GetkwError.h)

# List of sources
list(APPEND _sources_list Getkw.cpp GetkwError.cpp Section.cpp messages.cpp)

add_library(getkw OBJECT ${_sources_list} ${_headers_list})
if(BUILD_CUSTOM_BOOST)
  add_dependencies(getkw custom_boost)
endif()
set_target_properties(getkw PROPERTIES POSITION_INDEPENDENT_CODE 1
                                       INTERPROCEDURAL_OPTIMIZATION 1)

# Sets install directory for all the headers in the list
foreach(_header ${_headers_list})
   install(FILES ${_header} DESTINATION include/utils/getkw)
endforeach()
