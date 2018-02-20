add_library(utils-tests OBJECT
  utils_dipolar-potential.cpp
  utils_newton-potential.cpp
  )
set_target_properties(utils-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(utils-tests custom_boost)
endif()

# utils_newton-potential.cpp test
add_Catch_test(utils_newton-potential "utils;newton_potential;utils_newton-potential")

# utils_dipolar-potential.cpp test
add_Catch_test(utils_dipolar-potential "utils;dipolar_potential;utils_dipolar-potential")
