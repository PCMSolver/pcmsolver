configure_file(numerical_quadrature.cpp.in ${CMAKE_CURRENT_LIST_DIR}/numerical_quadrature.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/numerical_quadrature.cpp)
add_test(NAME numerical_quadrature COMMAND unit_tests.x [numerical_quadrature])

add_reference(molecule.npy ${TESTS_ROOT_DIRECTORY})
add_reference(molecule_1r.npy ${TESTS_ROOT_DIRECTORY})
