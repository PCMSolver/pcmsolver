# gepol_point_symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_point_symmetry.cpp)
add_Catch_test(gepol_point_symmetry "gepol;gepol_point_symmetry")

# gepol_point_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_point_from-file.cpp)
add_Catch_test(gepol_point_from-file "gepol;gepol_point_from-file")
set_tests_properties(gepol_point_from-file PROPERTIES DEPENDS gepol_point)

# gepol_CH3+_Cs.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_CH3+_Cs.cpp)
add_Catch_test(gepol_CH3+_Cs "gepol;gepol_CH3+_Cs")

# gepol_H3+_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_H3+_from-file.cpp)
add_Catch_test(gepol_H3+_from-file "gepol;gepol_H3+_from-file")
set_tests_properties(gepol_H3+_from-file PROPERTIES DEPENDS gepol_H3+)

# gepol_C2H4_D2h.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_C2H4_D2h.cpp)
add_Catch_test(gepol_C2H4_D2h "gepol;gepol_C2H4_D2h")

# gepol_H3+_C2v.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_H3+_C2v.cpp)
add_Catch_test(gepol_H3+_C2v "gepol;gepol_H3+_C2v")

# gepol_C6H6.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_C6H6.cpp)
add_Catch_test(gepol_C6H6 "gepol;gepol_C6H6")

# gepol_CO2_symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_CO2_symmetry.cpp)
add_Catch_test(gepol_CO2_symmetry "gepol;gepol_CO2_symmetry")

# gepol_NH3.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_NH3.cpp)
add_Catch_test(gepol_NH3 "gepol;gepol_NH3")

# gepol_point.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_point.cpp)
add_Catch_test(gepol_point "gepol;gepol_point")

# gepol_H3+.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_H3+.cpp)
add_Catch_test(gepol_H3+ "gepol;gepol_H3+")

# gepol_NH3_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/gepol_NH3_from-file.cpp)
add_Catch_test(gepol_NH3_from-file "gepol;gepol_NH3_from-file")
set_tests_properties(gepol_NH3_from-file PROPERTIES DEPENDS gepol_NH3)

