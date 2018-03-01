add_library(dielectric_profile-tests OBJECT
  one_layer.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(dielectric_profile-tests custom_boost)
endif()

# one_layer.cpp test
add_Catch_test(
  NAME
    one_layer
  LABELS
    dielectric_profile
    one_layer
 )

