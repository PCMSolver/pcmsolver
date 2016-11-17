# List of headers
list(APPEND _headers_list Anisotropic.hpp MembraneTanh.hpp Metal.hpp OneLayerErf.hpp OneLayerTanh.hpp ProfileForward.hpp ProfileTypes.hpp Sharp.hpp Uniform.hpp Yukawa.hpp)

# List of sources
list(APPEND _sources_list )

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
# Sets install directory for all the headers in the list
foreach(_header ${_headers_list})
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green/dielectric_profile)
endforeach()
