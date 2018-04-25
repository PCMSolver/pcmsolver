target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_anisotropic-gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_anisotropic-symmetry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_tanh-spherical-diffuse-gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_erf-spherical-diffuse-gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_log-spherical-diffuse-gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_sharp-gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_gepol-C2H4_D2h.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_gepol-NH3.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_gepol-NH3_from-file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_gepol-point_from-file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iefpcm_symmetry.cpp
  )

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
