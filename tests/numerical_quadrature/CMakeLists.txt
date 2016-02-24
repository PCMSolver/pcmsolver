add_library(numerical_quadrature-tests OBJECT numerical_quadrature.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(numerical_quadrature-tests custom_boost)
endif()

# numerical_quadrature.cpp test
configure_file(numerical_quadrature.cpp.in ${CMAKE_CURRENT_LIST_DIR}/numerical_quadrature.cpp @ONLY)
add_Catch_test(numerical_quadrature "numerical_quadrature")

