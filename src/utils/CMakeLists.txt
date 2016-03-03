# List of headers
list(APPEND headers_list Atom.hpp Factory.hpp ForId.hpp FortranCUtils.hpp Logger.hpp LoggerImpl.hpp MathUtils.hpp Molecule.hpp
                         PhysicalConstants.hpp QuadratureRules.hpp Solvent.hpp Sphere.hpp SplineFunction.hpp
                         Stencils.hpp Symmetry.hpp TestingMolecules.hpp Timer.hpp cnpy.hpp)
set_source_files_properties(cnpy.hpp PROPERTIES COMPILE_FLAGS -w)

# List of sources
list(APPEND sources_list Atom.cpp FortranCUtils.cpp Molecule.cpp Solvent.cpp Sphere.cpp Symmetry.cpp cnpy.cpp)
set_source_files_properties(cnpy.cpp PROPERTIES COMPILE_FLAGS -w)

add_library(utils OBJECT ${sources_list} ${headers_list})
set_target_properties(utils PROPERTIES POSITION_INDEPENDENT_CODE 1 INTERPROCEDURAL_OPTIMIZATION 1)
add_dependencies(utils generate-config-hpp)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/utils)
endforeach()

add_subdirectory(getkw)
