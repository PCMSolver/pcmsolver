add_library(green-tests OBJECT green_anisotropic_liquid.cpp green_ionic_liquid.cpp green_spherical_diffuse.cpp green_uniform_dielectric.cpp green_vacuum.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(green-tests custom_boost)
endif()

# green_spherical_diffuse.cpp test
add_Catch_test(green_spherical_diffuse "green;green_spherical_diffuse")

# green_vacuum.cpp test
add_Catch_test(green_vacuum "green;green_vacuum")

# green_ionic_liquid.cpp test
add_Catch_test(green_ionic_liquid "green;green_ionic_liquid")

# green_uniform_dielectric.cpp test
add_Catch_test(green_uniform_dielectric "green;green_uniform_dielectric")

# green_anisotropic_liquid.cpp test
add_Catch_test(green_anisotropic_liquid "green;green_anisotropic_liquid")
