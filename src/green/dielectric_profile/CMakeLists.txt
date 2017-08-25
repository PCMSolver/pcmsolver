# List of headers
list(APPEND _headers_list Anisotropic.hpp MembraneTanh.hpp Metal.hpp OneLayerErf.hpp OneLayerTanh.hpp ProfileTypes.hpp Sharp.hpp Uniform.hpp Yukawa.hpp)

# Sets install directory for all the headers in the list
foreach(_header ${_headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green/dielectric_profile)
endforeach()
