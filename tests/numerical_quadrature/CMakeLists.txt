target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/numerical_quadrature.cpp
  )

# numerical_quadrature.cpp test
add_Catch_test(
  NAME
    numerical_quadrature
  LABELS
    numerical_quadrature
  REFERENCE_FILES
    molecule_1r.npy
    molecule.npy
  )

