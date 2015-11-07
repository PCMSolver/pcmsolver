# green_spherical_diffuse.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/green_spherical_diffuse.cpp)
add_Catch_test(green_spherical_diffuse "green;green_spherical_diffuse")

# green_vacuum.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/green_vacuum.cpp)
add_Catch_test(green_vacuum "green;green_vacuum")

# green_ionic_liquid.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/green_ionic_liquid.cpp)
add_Catch_test(green_ionic_liquid "green;green_ionic_liquid")

# green_uniform_dielectric.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/green_uniform_dielectric.cpp)
add_Catch_test(green_uniform_dielectric "green;green_uniform_dielectric")

# green_anisotropic_liquid.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/green_anisotropic_liquid.cpp)
add_Catch_test(green_anisotropic_liquid "green;green_anisotropic_liquid")

