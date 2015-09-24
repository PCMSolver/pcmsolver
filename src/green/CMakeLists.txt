# List of headers
list(APPEND headers_list AnisotropicLiquid.hpp DerivativeTypes.hpp DerivativeUtils.hpp GreenUtils.hpp GreensFunction.hpp IGreensFunction.hpp InterfacesImpl.hpp IonicLiquid.hpp RegisterGreensFunctionToFactory.hpp SphericalDiffuse.hpp SphericalSharp.hpp UniformDielectric.hpp Vacuum.hpp)

# List of sources
list(APPEND sources_list )

# Update bar chart
if(BUILD_DOCUMENTATION)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR})
endif()

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/green)
endforeach()

add_subdirectory(dielectric_profile)
