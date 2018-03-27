add_library(cpcm-tests OBJECT
  cpcm_gepol-C2H4_D2h.cpp
  cpcm_gepol-H2O.cpp
  cpcm_gepol-NH3.cpp
  cpcm_gepol-NH3_from-file.cpp
  cpcm_gepol-point.cpp
  cpcm_gepol-point_from-file.cpp
  cpcm_symmetry.cpp
  cpcm_tsless-point.cpp
  cpcm_tsless-NH3.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(cpcm-tests custom_boost)
endif()

# cpcm_gepol-NH3_from-file.cpp test
add_Catch_test(
  NAME
    cpcm_gepol-NH3_from-file
  LABELS
    solver
    cpcm
    cpcm_gepol-NH3_from-file
  REFERENCE_FILES
    ASC-cpcm_gepol-NH3.npy
  DEPENDS
    cpcm_gepol-NH3
  )

# cpcm_gepol-H2O.cpp test
add_Catch_test(
  NAME
    cpcm_gepol-H2O
  LABELS
    solver
    cpcm
    cpcm_gepol-H2O
  REFERENCE_FILES
    ASC-cpcm_gepol-H2O.npy
  )

# cpcm_gepol-point_from-file.cpp test
add_Catch_test(
  NAME
    cpcm_gepol-point_from-file
  LABELS
    solver
    cpcm
    cpcm_gepol-point_from-file
  REFERENCE_FILES
    ASC-cpcm_gepol-point.npy
  DEPENDS
    cpcm_gepol-point
  )

# cpcm_symmetry.cpp test
add_Catch_test(
  NAME
    cpcm_symmetry
  LABELS
    solver
    cpcm
    cpcm_symmetry
  )

# cpcm_gepol-NH3.cpp test
add_Catch_test(
  NAME
    cpcm_gepol-NH3
  LABELS
    solver
    cpcm
    cpcm_gepol-NH3
  REFERENCE_FILES
    ASC-cpcm_gepol-NH3.npy
  )

# cpcm_gepol-point.cpp test
add_Catch_test(
  NAME
    cpcm_gepol-point
  LABELS
    solver
    cpcm
    cpcm_gepol-point
  REFERENCE_FILES
    ASC-cpcm_gepol-point.npy
  )

# cpcm_gepol-C2H4_D2h.cpp test
add_Catch_test(
  NAME
    cpcm_gepol-C2H4_D2h
  LABELS
    cpcm
    cpcm_symmetry
    cpcm_gepol-C2H4_D2h
  REFERENCE_FILES
    ASC-cpcm_gepol-C2H4_D2h.npy
  )

# FIXME Reference files for these two tests
# cpcm_tsless-point.cpp test
add_Catch_test(
  NAME
    cpcm_tsless-point
  LABELS
    cpcm
    cpcm_tsless-point
  )

# cpcm_tsless-NH3.cpp test
add_Catch_test(
  NAME
    cpcm_tsless-NH3
  LABELS
    cpcm
    cpcm_tsless-NH3
  )
