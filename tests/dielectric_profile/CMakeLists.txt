add_library(dielectric_profile-tests OBJECT
  one_layer.cpp
  )
set_target_properties(dielectric_profile-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(dielectric_profile-tests custom_boost)
endif()

# one_layer.cpp test
add_Catch_test(one_layer "dielectric_profile;one_layer")

