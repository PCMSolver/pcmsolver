add_library(utils-tests OBJECT utils_newton-potential.cpp utils_dipolar-potential.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(utils-tests custom_boost)
endif()

# utils_newton-potential.cpp test
add_Catch_test(utils_newton-potential "utils;newton_potential;utils_newton-potential")

# utils_dipolar-potential.cpp test
add_Catch_test(utils_dipolar-potential "utils;dipolar_potential;utils_dipolar-potential")
