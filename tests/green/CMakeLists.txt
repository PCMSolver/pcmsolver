add_library(green-tests OBJECT
  green_anisotropic_liquid.cpp
  green_ionic_liquid.cpp
  green_spherical_diffuse.cpp
  green_uniform_dielectric.cpp
  green_vacuum.cpp
  )
set_target_properties(green-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(green-tests custom_boost)
endif()

# green_spherical_diffuse.cpp test
add_Catch_test(green_spherical_diffuse "green;green_spherical_diffuse")
set_tests_properties(green_spherical_diffuse PROPERTIES COST 100.0)

# green_vacuum.cpp test
add_Catch_test(green_vacuum "green;green_vacuum")

# green_ionic_liquid.cpp test
add_Catch_test(green_ionic_liquid "green;green_ionic_liquid")

# green_uniform_dielectric.cpp test
add_Catch_test(green_uniform_dielectric "green;green_uniform_dielectric")

# green_anisotropic_liquid.cpp test
add_Catch_test(green_anisotropic_liquid "green;green_anisotropic_liquid")
