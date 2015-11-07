# input.cpp test
configure_file(input.cpp.in ${CMAKE_CURRENT_LIST_DIR}/input.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/input.cpp)
add_Catch_test(input "input")

