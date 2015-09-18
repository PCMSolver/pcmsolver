set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/surface_function.cpp)
add_test(NAME surface_function COMMAND unit_tests.x surface_function)
