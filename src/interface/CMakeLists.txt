# List of headers
list(APPEND headers_list Meddle.hpp Input.hpp)

# List of sources
list(APPEND sources_list Meddle.cpp Input.cpp)

add_library(interface OBJECT ${sources_list} ${headers_list})
set_target_properties(interface PROPERTIES POSITION_INDEPENDENT_CODE 1 INTERPROCEDURAL_OPTIMIZATION 1)
add_dependencies(interface generate-config-hpp)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION include/interface)
endforeach()
