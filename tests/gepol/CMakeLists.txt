add_library(gepol-tests OBJECT
  gepol_C2H4_D2h.cpp
  gepol_C6H6.cpp
  gepol_CH3+_Cs.cpp
  gepol_CO2_symmetry.cpp
  gepol_H3+.cpp
  gepol_H3+_C2v.cpp
  gepol_H3+_from-file.cpp
  gepol_NH3.cpp
  gepol_NH3_from-file.cpp
  gepol_point.cpp
  gepol_point_from-file.cpp
  gepol_point_symmetry.cpp
  )
set_target_properties(gepol-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(gepol-tests custom_boost)
endif()

# gepol_point_symmetry.cpp test
add_Catch_test(gepol_point_symmetry "gepol;gepol_point_symmetry")

# gepol_point_from-file.cpp test
add_Catch_test(gepol_point_from-file "gepol;gepol_point_from-file")
set_tests_properties(gepol_point_from-file PROPERTIES DEPENDS gepol_point)

# gepol_CH3+_Cs.cpp test
add_Catch_test(gepol_CH3+_Cs "gepol;gepol_CH3+_Cs")

# gepol_H3+_from-file.cpp test
add_Catch_test(gepol_H3+_from-file "gepol;gepol_H3+_from-file")
set_tests_properties(gepol_H3+_from-file PROPERTIES DEPENDS gepol_H3+)

# gepol_C2H4_D2h.cpp test
add_Catch_test(gepol_C2H4_D2h "gepol;gepol_C2H4_D2h")

# gepol_H3+_C2v.cpp test
add_Catch_test(gepol_H3+_C2v "gepol;gepol_H3+_C2v")

# gepol_C6H6.cpp test
add_Catch_test(gepol_C6H6 "gepol;gepol_C6H6")

# gepol_CO2_symmetry.cpp test
add_Catch_test(gepol_CO2_symmetry "gepol;gepol_CO2_symmetry")

# gepol_NH3.cpp test
add_Catch_test(gepol_NH3 "gepol;gepol_NH3")

# gepol_point.cpp test
add_Catch_test(gepol_point "gepol;gepol_point")

# gepol_H3+.cpp test
add_Catch_test(gepol_H3+ "gepol;gepol_H3+")

# gepol_NH3_from-file.cpp test
add_Catch_test(gepol_NH3_from-file "gepol;gepol_NH3_from-file")
set_tests_properties(gepol_NH3_from-file PROPERTIES DEPENDS gepol_NH3)
