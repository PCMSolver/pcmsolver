# List of headers
list(APPEND headers_list Atom.hpp CavityData.hpp Citation.hpp Cxx11Workarounds.hpp ErrorHandling.hpp Factory.hpp ForId.hpp FortranCUtils.hpp GreenData.hpp Input.hpp Interpolation.hpp Logger.hpp LoggerImpl.hpp LoggerInterface.hpp MathUtils.hpp Molecule.hpp PhysicalConstants.hpp QuadratureRules.hpp RadialFunction.hpp Solvent.hpp SolverData.hpp Sphere.hpp SplineFunction.hpp Stencils.hpp SurfaceFunction.hpp Symmetry.hpp TestingMolecules.hpp Timer.hpp TimerInterface.hpp Vector2.hpp Vector3.hpp cnpy.hpp)
set_source_files_properties(cnpy.hpp PROPERTIES COMPILE_FLAGS -w)

# List of sources
list(APPEND sources_list Atom.cpp FortranCUtils.cpp Input.cpp Interpolation.cpp Molecule.cpp Solvent.cpp Sphere.cpp SurfaceFunction.cpp Symmetry.cpp cnpy.cpp)
set_source_files_properties(cnpy.cpp PROPERTIES COMPILE_FLAGS -w)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/utils)
endforeach()
