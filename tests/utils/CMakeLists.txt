target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/utils_dipolar-potential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils_newton-potential.cpp
  )

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
