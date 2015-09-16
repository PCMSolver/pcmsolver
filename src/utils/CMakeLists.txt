# List of headers
list(APPEND headers_list Atom.hpp CavityData.hpp Citation.hpp Cxx11Workarounds.hpp ErrorHandling.hpp Exception.hpp Factory.hpp ForId.hpp GeneralPurpose.hpp GreenData.hpp Input.hpp Interpolation.hpp Logger.hpp LoggerImpl.hpp LoggerInterface.hpp MathUtils.hpp Molecule.hpp PhysicalConstants.hpp QuadratureRules.hpp RadialFunction.hpp Solvent.hpp SolverData.hpp Sphere.hpp SplineFunction.hpp Stencils.hpp SurfaceFunction.hpp Symmetry.hpp TestingMolecules.hpp Timer.hpp TimerInterface.hpp Vector2.hpp Vector3.hpp cnpy.hpp cnpyPimpl.hpp)

# List of sources
list(APPEND sources_list Atom.cpp Exception.cpp Input.cpp Interpolation.cpp Molecule.cpp Solvent.cpp Sphere.cpp SurfaceFunction.cpp Symmetry.cpp cnpy.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(utils STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(utils libgetkw)
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES utils)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(utils custom_boost)
endif()
install(TARGETS utils ARCHIVE DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/utils)
endforeach()
