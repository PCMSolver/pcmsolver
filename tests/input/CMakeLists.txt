add_library(input-tests OBJECT input.cpp)

# input.cpp test
configure_file(input.cpp.in ${CMAKE_CURRENT_LIST_DIR}/input.cpp @ONLY)
add_Catch_test(input "input")

