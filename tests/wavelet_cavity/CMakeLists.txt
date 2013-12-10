get_filename_component(test_name wavelet_cavity.cpp NAME_WE)
add_executable(${test_name}.x wavelet_cavity.cpp)      	
add_dependencies(${test_name}.x googletest)
target_link_libraries(${test_name}.x
	                  solver green cavity wem pwl wavcav utils
                      ${GTEST_LIBS_DIR}/libgtest.a
                      ${GTEST_LIBS_DIR}/libgtest_main.a
                      ${LIBS}
                     )
add_test(NAME ${test_name} COMMAND ${test_name}.x)
