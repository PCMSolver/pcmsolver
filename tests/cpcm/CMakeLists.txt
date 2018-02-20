add_library(cpcm-tests OBJECT
  cpcm_gepol-C2H4_D2h.cpp
  cpcm_gepol-H2O.cpp
  cpcm_gepol-NH3.cpp
  cpcm_gepol-NH3_from-file.cpp
  cpcm_gepol-point.cpp
  cpcm_gepol-point_from-file.cpp
  cpcm_symmetry.cpp
  )
set_target_properties(cpcm-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(cpcm-tests custom_boost)
endif()

# Copy reference files to ${PROJECT_BINARY_DIR}/tests/bi_operators (aka ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND reference_files
  ASC-cpcm_gepol-point.npy
  ASC-cpcm_gepol-NH3.npy
  ASC-cpcm_gepol-H2O.npy
  ASC-cpcm_gepol-C2H4_D2h.npy
  )
file(COPY ${reference_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# cpcm_gepol-NH3_from-file.cpp test
add_Catch_test(cpcm_gepol-NH3_from-file "solver;cpcm;cpcm_gepol-NH3_from-file")
set_tests_properties(cpcm_gepol-NH3_from-file PROPERTIES DEPENDS cpcm_gepol-NH3)

# cpcm_gepol-H2O.cpp test
add_Catch_test(cpcm_gepol-H2O "solver;cpcm;cpcm_gepol-H2O")

# cpcm_gepol-point_from-file.cpp test
add_Catch_test(cpcm_gepol-point_from-file "solver;cpcm;cpcm_gepol-point_from-file")
set_tests_properties(cpcm_gepol-point_from-file PROPERTIES DEPENDS cpcm_gepol-point)

# cpcm_symmetry.cpp test
add_Catch_test(cpcm_symmetry "solver;cpcm;cpcm_symmetry")

# cpcm_gepol-NH3.cpp test
add_Catch_test(cpcm_gepol-NH3 "solver;cpcm;cpcm_gepol-NH3")

# cpcm_gepol-point.cpp test
add_Catch_test(cpcm_gepol-point "solver;cpcm;cpcm_gepol-point")

# cpcm_gepol-C2H4_D2h.cpp test
add_Catch_test(cpcm_gepol-C2H4_D2h "cpcm;cpcm_symmetry;cpcm_gepol-C2H4_D2h")
