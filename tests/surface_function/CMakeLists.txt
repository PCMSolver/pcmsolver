add_library(surface_function-tests OBJECT surface_function.cpp)

# surface_function.cpp test
add_Catch_test(surface_function "surface_function")

