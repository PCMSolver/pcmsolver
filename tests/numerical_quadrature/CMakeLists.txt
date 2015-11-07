# numerical_quadrature.cpp test
configure_file(numerical_quadrature.cpp.in ${CMAKE_CURRENT_LIST_DIR}/numerical_quadrature.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/numerical_quadrature.cpp)
add_Catch_test(numerical_quadrature "numerical_quadrature")

