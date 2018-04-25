# List of headers
list(APPEND headers_list
  Anisotropic.hpp
  Metal.hpp
  OneLayerErf.hpp
  OneLayerLog.hpp
  OneLayerTanh.hpp
  Sharp.hpp
  Uniform.hpp
  Yukawa.hpp
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS _headers_list)
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green/dielectric_profile)
endforeach()
