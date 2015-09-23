configure_file(bi_operators_collocation.cpp.in ${CMAKE_CURRENT_LIST_DIR}/bi_operators_collocation.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/bi_operators_collocation.cpp)
add_test(NAME bi_operators_collocation COMMAND unit_tests.x [bi_operators_collocation])

configure_file(bi_operators_numerical.cpp.in ${CMAKE_CURRENT_LIST_DIR}/bi_operators_numerical.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/bi_operators_numerical.cpp)
add_test(NAME bi_operators_numerical COMMAND unit_tests.x [bi_operators_numerical])

configure_file(bi_operators_purisima.cpp.in ${CMAKE_CURRENT_LIST_DIR}/bi_operators_purisima.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/bi_operators_purisima.cpp)
add_test(NAME bi_operators_purisima COMMAND unit_tests.x [bi_operators_purisima])

# This executable updates the .npy files containing the reference values
# The reference are to be copied into the source directory for the tests
list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
add_executable(update_reference_files.x update_reference_files.cpp)
target_link_libraries(update_reference_files.x pcm ${external_libraries})

