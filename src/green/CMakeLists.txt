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

# List of sources
list(APPEND sources_list
AnisotropicLiquid.cpp
IonicLiquid.cpp
SphericalDiffuse.cpp
UniformDielectric.cpp
Vacuum.cpp
)

add_library(green OBJECT ${sources_list} ${headers_list})
target_compile_definitions(green
  PUBLIC
    PCMSolver_EXPORTS
  )
target_compile_options(green
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )
target_include_directories(green
  PUBLIC
    ${PROJECT_SOURCE_DIR}/src/dielectric_profile
  )
set_target_properties(green
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(green custom_boost)
endif()
add_dependencies(green generate-config-hpp)

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green)
endforeach()

add_subdirectory(dielectric_profile)
