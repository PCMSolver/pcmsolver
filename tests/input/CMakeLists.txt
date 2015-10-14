configure_file(input.cpp.in ${CMAKE_CURRENT_LIST_DIR}/input.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/input.cpp)
add_test(NAME input COMMAND unit_tests.x [input])
