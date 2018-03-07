add_library(numerical_quadrature-tests OBJECT
  numerical_quadrature.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(numerical_quadrature-tests custom_boost)
endif()

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

