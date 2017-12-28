# List of headers
list(APPEND headers_list
Input.hpp
Meddle.hpp
)

# List of sources
list(APPEND sources_list
Input.cpp
Meddle.cpp
PhysicalConstants.cpp
)

add_library(interface OBJECT ${sources_list} ${headers_list})
set_target_properties(interface PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(interface generate-config-hpp)
target_compile_options(interface PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/interface)
endforeach()
