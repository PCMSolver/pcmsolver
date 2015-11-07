# one_layer.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/one_layer.cpp)
add_Catch_test(one_layer "dielectric_profile;one_layer")

