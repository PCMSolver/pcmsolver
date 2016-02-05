add_library(dielectric_profile-tests OBJECT one_layer.cpp)

# one_layer.cpp test
add_Catch_test(one_layer "dielectric_profile;one_layer")

