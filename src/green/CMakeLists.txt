target_sources(pcm-objlib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/AnisotropicLiquid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IonicLiquid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SphericalDiffuse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UniformDielectric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Vacuum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SphericalSharp.cpp
  )

# List of headers
list(APPEND headers_list
  AnisotropicLiquid.hpp
  DerivativeTypes.hpp
  Green.hpp
  GreenData.hpp
  GreenUtils.hpp
  GreensFunction.hpp
  IGreensFunction.hpp
  InterfacesImpl.hpp
  IonicLiquid.hpp
  SphericalDiffuse.hpp
  SphericalSharp.hpp
  UniformDielectric.hpp
  Vacuum.hpp
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS headers_list)
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green)
endforeach()

add_subdirectory(dielectric_profile)
