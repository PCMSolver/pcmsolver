add_library(iefpcm-tests OBJECT
  iefpcm_anisotropic-gepol-point.cpp
  iefpcm_anisotropic-symmetry.cpp
  iefpcm_tanh-spherical-diffuse-gepol-point.cpp
  iefpcm_erf-spherical-diffuse-gepol-point.cpp
  iefpcm_log-spherical-diffuse-gepol-point.cpp
  iefpcm_sharp-gepol-point.cpp
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
add_Catch_test(
  NAME
    iefpcm_gepol-NH3
  LABELS
    solver
    iefpcm
    iefpcm_gepol-NH3
  )

# iefpcm_anisotropic-symmetry.cpp test
add_Catch_test(
  NAME
    iefpcm_anisotropic-symmetry
  LABELS
    solver
    iefpcm
    iefpcm_anisotropic-symmetry
    anisotropic
  )

# iefpcm_tanh-spherical-diffuse-gepol-point.cpp test
add_Catch_test(
  NAME
    iefpcm_tanh-spherical-diffuse-gepol-point
  LABELS
    solver
    iefpcm
    iefpcm_tanh-spherical-diffuse-gepol-point
    anisotropic
  COST
    200.0
  )

# iefpcm_erf-spherical-diffuse-gepol-point.cpp test
add_Catch_test(
  NAME
    iefpcm_erf-spherical-diffuse-gepol-point
  LABELS
    solver
    iefpcm
    iefpcm_erf-spherical-diffuse-gepol-point
    anisotropic
  COST
    200.0
  )

# iefpcm_log-spherical-diffuse-gepol-point.cpp test
add_Catch_test(
  NAME
    iefpcm_log-spherical-diffuse-gepol-point
  LABELS
    solver
    iefpcm
    iefpcm_log-spherical-diffuse-gepol-point
    anisotropic
  COST
    200.0
  )

# iefpcm_sharp-gepol-point.cpp test
add_Catch_test(
  NAME
    iefpcm_sharp-gepol-point
  LABELS
    solver
    iefpcm
    iefpcm_sharp-gepol-point
    anisotropic
  )

# iefpcm_gepol-point.cpp test
add_Catch_test(
  NAME
    iefpcm_gepol-point
  LABELS
    solver
    iefpcm
    iefpcm_gepol-point
  )

# iefpcm_symmetry.cpp test
add_Catch_test(
  NAME
    iefpcm_symmetry
  LABELS
    solver
    iefpcm
    iefpcm_symmetry
  )

# iefpcm_gepol-C2H4_D2h.cpp test
add_Catch_test(
  NAME
    iefpcm_gepol-C2H4_D2h
  LABELS
    iefpcm
    iefpcm_symmetry
    iefpcm_gepol-C2H4_D2h
  )

# iefpcm_gepol-NH3_from-file.cpp test
add_Catch_test(
  NAME
    iefpcm_gepol-NH3_from-file
  LABELS
    solver
    iefpcm
    iefpcm_gepol-NH3_from-file
  DEPENDS
    iefpcm_gepol-NH3
  )

# iefpcm_anisotropic-gepol-point.cpp test
add_Catch_test(
  NAME
    iefpcm_anisotropic-gepol-point
  LABELS
    solver
    iefpcm
    iefpcm_anisotropic-gepol-point
    anisotropic
  )

# iefpcm_gepol-point_from-file.cpp test
add_Catch_test(
  NAME
    iefpcm_gepol-point_from-file
  LABELS
    solver
    iefpcm
    iefpcm_gepol-point_from-file
  DEPENDS
    iefpcm_gepol-point
  )
