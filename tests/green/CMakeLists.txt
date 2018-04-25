target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/green_anisotropic_liquid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/green_ionic_liquid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/green_spherical_diffuse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/green_uniform_dielectric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/green_vacuum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/green_spherical_sharp.cpp
  )

# green_spherical_diffuse.cpp test
add_Catch_test(
  NAME
    green_spherical_diffuse
  LABELS
    green
    green_spherical_diffuse
  COST
    100.0
  )

# green_vacuum.cpp test
add_Catch_test(
  NAME
    green_vacuum
  LABELS
    green
    green_vacuum
  )

# green_ionic_liquid.cpp test
add_Catch_test(
  NAME
    green_ionic_liquid
  LABELS
    green
    green_ionic_liquid
  )

# green_uniform_dielectric.cpp test
add_Catch_test(
  NAME
    green_uniform_dielectric
  LABELS
    green
    green_uniform_dielectric
  )

# green_anisotropic_liquid.cpp test
add_Catch_test(
  NAME
    green_anisotropic_liquid
  LABELS
    green
    green_anisotropic_liquid
  )

# green_sharp_diffuse.cpp test
add_Catch_test(
  NAME
    green_sharp_diffuse
  LABELS
    green
    green_sharp_diffuse
  COST
    100.0
  )
