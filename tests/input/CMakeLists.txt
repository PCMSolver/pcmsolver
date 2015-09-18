configure_file(input.cpp.in ${CMAKE_CURRENT_LIST_DIR}/input.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/input.cpp)
add_test(NAME input COMMAND unit_tests.x [input])

add_reference(@tsless.inp ${TESTS_ROOT_DIRECTORY})
add_reference(@restart.inp ${TESTS_ROOT_DIRECTORY})
add_reference(@wavelet.inp ${TESTS_ROOT_DIRECTORY})
add_reference(@diffuse.inp ${TESTS_ROOT_DIRECTORY})
