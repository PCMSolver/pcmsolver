add_library(numerical_quadrature-tests OBJECT
  numerical_quadrature.cpp
  )
set_target_properties(numerical_quadrature-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(numerical_quadrature-tests custom_boost)
endif()

# Copy reference files to ${PROJECT_BINARY_DIR}/tests/numerical_quadrature (aka ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND reference_files
  molecule_1r.npy
  molecule.npy
  )
file(COPY ${reference_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# numerical_quadrature.cpp test
add_Catch_test(numerical_quadrature "numerical_quadrature")

