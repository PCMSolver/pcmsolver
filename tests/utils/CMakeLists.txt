add_library(utils-tests OBJECT
  utils_dipolar-potential.cpp
  utils_newton-potential.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(utils-tests custom_boost)
endif()

# utils_newton-potential.cpp test
add_Catch_test(
  NAME
    utils_newton-potential
  LABELS
    utils
    newton_potential
    utils_newton-potential
  )

# utils_dipolar-potential.cpp test
add_Catch_test(
  NAME
    utils_dipolar-potential
  LABELS
    utils
    dipolar_potential
    utils_dipolar-potential
  )
