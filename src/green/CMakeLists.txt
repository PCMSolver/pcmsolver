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

include_directories(dielectric_profile)
add_library(green OBJECT ${sources_list} ${headers_list})
set_target_properties(green PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(green generate-config-hpp)
target_compile_options(green PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/green)
endforeach()

add_subdirectory(dielectric_profile)
