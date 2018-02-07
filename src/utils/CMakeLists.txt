# List of headers
list(APPEND headers_list
Atom.hpp
ChargeDistribution.hpp
Factory.hpp
Legendre.hpp
Logger.hpp
LoggerImpl.hpp
MathUtils.hpp
Molecule.hpp
QuadratureRules.hpp
Solvent.hpp
Sphere.hpp
SplineFunction.hpp
Stencils.hpp
Symmetry.hpp
Timer.hpp
cnpy.hpp
)

# List of sources
list(APPEND sources_list
Atom.cpp
ChargeDistribution.cpp
Molecule.cpp
QuadratureRules.cpp
Solvent.cpp
Sphere.cpp
Symmetry.cpp
cnpy.cpp
)

add_library(utils OBJECT ${sources_list} ${headers_list})
set_target_properties(utils
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
add_dependencies(utils generate-config-hpp)
target_compile_options(utils
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils)
endforeach()

add_subdirectory(getkw)
