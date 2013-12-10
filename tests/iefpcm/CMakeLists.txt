get_filename_component(test_name iefpcm-gepol.cpp NAME_WE)
add_executable(${test_name}.x iefpcm-gepol.cpp)      	
add_dependencies(${test_name}.x googletest)
target_link_libraries(${test_name}.x
	                  solver green cavity pedra utils
                      ${GTEST_LIBS_DIR}/libgtest.a
                      ${GTEST_LIBS_DIR}/libgtest_main.a
                      ${LIBS}
                     )
add_test(NAME ${test_name} COMMAND ${test_name}.x)

get_filename_component(test_name iefpcm-tsless.cpp NAME_WE)
add_executable(${test_name}.x iefpcm-tsless.cpp)      	
add_dependencies(${test_name}.x googletest)
target_link_libraries(${test_name}.x
	                  solver green cavity tsless pedra utils
                      ${GTEST_LIBS_DIR}/libgtest.a
                      ${GTEST_LIBS_DIR}/libgtest_main.a
                      ${LIBS}
                     )
add_test(NAME ${test_name} COMMAND ${test_name}.x)
