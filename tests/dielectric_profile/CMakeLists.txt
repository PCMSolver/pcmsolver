target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/one_layer.cpp
  )

# one_layer.cpp test
add_Catch_test(
  NAME
    one_layer
  LABELS
    dielectric_profile
    one_layer
 )
