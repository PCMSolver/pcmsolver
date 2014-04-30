# This macro is used to add a unit test using Google Unit Testing framework
macro(add_googletest test_name my_libraries external_libraries)
    get_filename_component(the_name ${test_name} NAME_WE)
    add_executable(${the_name}.x ${the_name}.cpp)      	
    add_dependencies(${the_name}.x googletest)
    target_link_libraries(${the_name}.x
                          ${my_libraries}
                          ${GTEST_LIBS_DIR}/libgtest.a
                          ${GTEST_LIBS_DIR}/libgtest_main.a
                          ${external_libraries}
                         )
    add_test(NAME ${the_name} COMMAND ${the_name}.x)
endmacro()

# This macro is used to add a unit test using Boost Unit Testing framework
macro(add_boosttest test_name my_libraries external_libraries)
    get_filename_component(the_name ${test_name} NAME_WE)
    add_executable(${the_name}.x ${the_name}.cpp)      	
    target_link_libraries(${the_name}.x
                          ${my_libraries}
			  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                          ${external_libraries}
                         )
    add_test(NAME ${the_name} COMMAND ${the_name}.x)
endmacro()
