target_sources(pcm-objlib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Atom.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ChargeDistribution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QuadratureRules.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Solvent.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Sphere.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Symmetry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cnpy.cpp
  )

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
# Sets install directory for all the headers in the list
foreach(_header IN LISTS headers_list)
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils)
endforeach()

add_subdirectory(getkw)
