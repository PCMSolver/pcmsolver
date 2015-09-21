set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/one_layer.cpp)
add_test(NAME one_layer COMMAND unit_tests.x [one_layer])
