add_library(tsless-tests OBJECT
  tsless_HF.cpp
  tsless_NH3.cpp
  tsless_point.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(tsless-tests custom_boost)
endif()

# tsless_HF.cpp test
add_Catch_test(
  NAME
    tsless_HF
  LABELS
    tsless
    tsless_HF
  )

# tsless_point.cpp test
add_Catch_test(
  NAME
    tsless_point
  LABELS
    tsless
    tsless_point
  )

# tsless_NH3.cpp test
add_Catch_test(
  NAME
    tsless_NH3
  LABELS
    tsless
    tsless_NH3
  )

