# List of headers
list(APPEND _headers_list Anisotropic.hpp MembraneTanh.hpp Metal.hpp OneLayerErf.hpp OneLayerTanh.hpp ProfileTypes.hpp Sharp.hpp Uniform.hpp Yukawa.hpp)

# List of sources
list(APPEND _sources_list )

add_library(dielectric_profile OBJECT ${_headers_list} ${_sources_list})
set_target_properties(dielectric_profile PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(dielectric_profile generate-config-hpp)
target_compile_options(dielectric_profile PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")
# Sets install directory for all the headers in the list
foreach(_header ${_headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green/dielectric_profile)
endforeach()
