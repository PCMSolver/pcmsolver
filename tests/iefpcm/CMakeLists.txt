# iefpcm_gepol-NH3.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_gepol-NH3.cpp)
add_Catch_test(iefpcm_gepol-NH3 "solver;iefpcm;iefpcm_gepol-NH3")

# iefpcm_anisotropic-symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_anisotropic-symmetry.cpp)
add_Catch_test(iefpcm_anisotropic-symmetry "solver;iefpcm;iefpcm_anisotropic-symmetry;anisotropic")

# iefpcm_diffuse-gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_diffuse-gepol-point.cpp)
add_Catch_test(iefpcm_diffuse-gepol-point "solver;iefpcm;iefpcm_diffuse-gepol-point;anisotropic")

# iefpcm_gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_gepol-point.cpp)
add_Catch_test(iefpcm_gepol-point "solver;iefpcm;iefpcm_gepol-point")

# iefpcm_symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_symmetry.cpp)
add_Catch_test(iefpcm_symmetry "solver;iefpcm;iefpcm_symmetry")

# iefpcm_gepol-C2H4_D2h.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_gepol-C2H4_D2h.cpp)
add_Catch_test(iefpcm_gepol-C2H4_D2h "iefpcm;iefpcm_symmetry;iefpcm_gepol-C2H4_D2h")

# iefpcm_gepol-NH3_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_gepol-NH3_from-file.cpp)
add_Catch_test(iefpcm_gepol-NH3_from-file "solver;iefpcm;iefpcm_gepol-NH3_from-file")
set_tests_properties(iefpcm_gepol-NH3_from-file PROPERTIES DEPENDS iefpcm_gepol-NH3)

# iefpcm_anisotropic-gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_anisotropic-gepol-point.cpp)
add_Catch_test(iefpcm_anisotropic-gepol-point "solver;iefpcm;iefpcm_anisotropic-gepol-point;anisotropic")

# iefpcm_gepol-point_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/iefpcm_gepol-point_from-file.cpp)
add_Catch_test(iefpcm_gepol-point_from-file "solver;iefpcm;iefpcm_gepol-point_from-file")
set_tests_properties(iefpcm_gepol-point_from-file PROPERTIES DEPENDS iefpcm_gepol-point)

