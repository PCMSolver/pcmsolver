# cpcm_gepol-NH3_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/cpcm_gepol-NH3_from-file.cpp)
add_Catch_test(cpcm_gepol-NH3_from-file "solver;cpcm;cpcm_gepol-NH3_from-file")
set_tests_properties(cpcm_gepol-NH3_from-file PROPERTIES DEPENDS cpcm_gepol-NH3)

# cpcm_gepol-point_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/cpcm_gepol-point_from-file.cpp)
add_Catch_test(cpcm_gepol-point_from-file "solver;cpcm;cpcm_gepol-point_from-file")
set_tests_properties(cpcm_gepol-point_from-file PROPERTIES DEPENDS cpcm_gepol-point)

# cpcm_symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/cpcm_symmetry.cpp)
add_Catch_test(cpcm_symmetry "solver;cpcm;cpcm_symmetry")

# cpcm_gepol-NH3.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/cpcm_gepol-NH3.cpp)
add_Catch_test(cpcm_gepol-NH3 "solver;cpcm;cpcm_gepol-NH3")

# cpcm_gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/cpcm_gepol-point.cpp)
add_Catch_test(cpcm_gepol-point "solver;cpcm;cpcm_gepol-point")

# cpcm_gepol-C2H4_D2h.cpp test
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/cpcm_gepol-C2H4_D2h.cpp)
add_Catch_test(cpcm_gepol-C2H4_D2h "cpcm;cpcm_symmetry;cpcm_gepol-C2H4_D2h")

