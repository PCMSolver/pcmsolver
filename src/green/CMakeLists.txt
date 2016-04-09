# List of headers
list(APPEND headers_list AnisotropicLiquid.hpp DerivativeTypes.hpp DerivativeUtils.hpp GreenData.hpp GreenUtils.hpp GreensFunction.hpp
                         IGreensFunction.hpp InterfacesImpl.hpp IonicLiquid.hpp
                         RegisterGreensFunctionToFactory.hpp SphericalDiffuse.hpp
                         SphericalSharp.hpp UniformDielectric.hpp Vacuum.hpp)

# List of sources
list(APPEND sources_list )

add_library(green OBJECT ${sources_list} ${headers_list})
set_target_properties(green PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(green generate-config-hpp)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/green)
endforeach()

add_subdirectory(dielectric_profile)
