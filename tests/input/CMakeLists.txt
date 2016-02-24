add_library(input-tests OBJECT input.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(input-tests custom_boost)
endif()

# input.cpp test
configure_file(input.cpp.in ${CMAKE_CURRENT_LIST_DIR}/input.cpp @ONLY)
add_Catch_test(input "input")
