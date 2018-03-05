add_library(iefpcm-tests OBJECT
  iefpcm_anisotropic-gepol-point.cpp
  iefpcm_anisotropic-symmetry.cpp
  iefpcm_tanh-spherical-diffuse-gepol-point.cpp
  iefpcm_erf-spherical-diffuse-gepol-point.cpp
  iefpcm_log-spherical-diffuse-gepol-point.cpp
  iefpcm_gepol-C2H4_D2h.cpp
  iefpcm_gepol-NH3.cpp
  iefpcm_gepol-NH3_from-file.cpp
  iefpcm_gepol-point.cpp
  iefpcm_gepol-point_from-file.cpp
  iefpcm_symmetry.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(iefpcm-tests custom_boost)
endif()

# iefpcm_gepol-NH3.cpp test
add_Catch_test(iefpcm_gepol-NH3 "solver;iefpcm;iefpcm_gepol-NH3")

# iefpcm_anisotropic-symmetry.cpp test
add_Catch_test(iefpcm_anisotropic-symmetry "solver;iefpcm;iefpcm_anisotropic-symmetry;anisotropic")

# iefpcm_tanh-spherical-diffuse-gepol-point.cpp test
add_Catch_test(iefpcm_tanh-spherical-diffuse-gepol-point
  "solver;iefpcm;iefpcm_tanh-spherical-diffuse-gepol-point;anisotropic")
set_tests_properties(iefpcm_tanh-spherical-diffuse-gepol-point PROPERTIES COST 200.0)

# iefpcm_erf-spherical-diffuse-gepol-point.cpp test
add_Catch_test(iefpcm_erf-spherical-diffuse-gepol-point
  "solver;iefpcm;iefpcm_erf-spherical-diffuse-gepol-point;anisotropic")
set_tests_properties(iefpcm_erf-spherical-diffuse-gepol-point PROPERTIES COST 200.0)

# iefpcm_log-spherical-diffuse-gepol-point.cpp test
add_Catch_test(iefpcm_log-spherical-diffuse-gepol-point
  "solver;iefpcm;iefpcm_log-spherical-diffuse-gepol-point;anisotropic")
set_tests_properties(iefpcm_log-spherical-diffuse-gepol-point PROPERTIES COST 200.0)

# iefpcm_gepol-point.cpp test
add_Catch_test(iefpcm_gepol-point "solver;iefpcm;iefpcm_gepol-point")

# iefpcm_symmetry.cpp test
add_Catch_test(iefpcm_symmetry "solver;iefpcm;iefpcm_symmetry")

# iefpcm_gepol-C2H4_D2h.cpp test
add_Catch_test(iefpcm_gepol-C2H4_D2h "iefpcm;iefpcm_symmetry;iefpcm_gepol-C2H4_D2h")

# iefpcm_gepol-NH3_from-file.cpp test
add_Catch_test(iefpcm_gepol-NH3_from-file "solver;iefpcm;iefpcm_gepol-NH3_from-file")
set_tests_properties(iefpcm_gepol-NH3_from-file PROPERTIES DEPENDS iefpcm_gepol-NH3)

# iefpcm_anisotropic-gepol-point.cpp test
add_Catch_test(iefpcm_anisotropic-gepol-point "solver;iefpcm;iefpcm_anisotropic-gepol-point;anisotropic")

# iefpcm_gepol-point_from-file.cpp test
add_Catch_test(iefpcm_gepol-point_from-file "solver;iefpcm;iefpcm_gepol-point_from-file")
set_tests_properties(iefpcm_gepol-point_from-file PROPERTIES DEPENDS iefpcm_gepol-point)

