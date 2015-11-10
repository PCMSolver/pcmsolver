# bi_operators_purisima.cpp test
configure_file(bi_operators_purisima.cpp.in ${CMAKE_CURRENT_LIST_DIR}/bi_operators_purisima.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/bi_operators_purisima.cpp)
add_Catch_test(bi_operators_purisima "bi_operators;bi_operators_purisima")

# bi_operators_numerical.cpp test
configure_file(bi_operators_numerical.cpp.in ${CMAKE_CURRENT_LIST_DIR}/bi_operators_numerical.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/bi_operators_numerical.cpp)
add_Catch_test(bi_operators_numerical "bi_operators;bi_operators_numerical")

# bi_operators_collocation.cpp test
configure_file(bi_operators_collocation.cpp.in ${CMAKE_CURRENT_LIST_DIR}/bi_operators_collocation.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/bi_operators_collocation.cpp)
add_Catch_test(bi_operators_collocation "bi_operators;bi_operators_collocation")

# This executable updates the .npy files containing the reference values
if(BUILD_STANDALONE)
  list(APPEND external_libraries ${ZLIB_LIBRARIES})
  add_executable(update_reference_files update_reference_files.cpp)
  target_link_libraries(update_reference_files pcm ${external_libraries})
endif()
