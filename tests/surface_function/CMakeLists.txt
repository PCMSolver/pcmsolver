add_library(surface_function-tests OBJECT surface_function.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(surface_function-tests custom_boost)
endif()

# surface_function.cpp test
add_Catch_test(surface_function "surface_function")

