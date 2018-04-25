target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_C2H4_D2h.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_C6H6.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_CH3+_Cs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_CO2_symmetry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_H3+.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_H3+_C2v.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_H3+_from-file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_NH3.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_NH3_from-file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_point_from-file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gepol_point_symmetry.cpp
  )

# gepol_point_symmetry.cpp test
add_Catch_test(
  NAME
    gepol_point_symmetry
  LABELS
    gepol
    gepol_point_symmetry
  )

# gepol_point_from-file.cpp test
add_Catch_test(
  NAME
    gepol_point_from-file
  LABELS
    gepol
    gepol_point_from-file
  DEPENDS
    gepol_point
  )

# gepol_CH3+_Cs.cpp test
add_Catch_test(
  NAME
    gepol_CH3+_Cs
  LABELS
    gepol
    gepol_CH3+_Cs
  )

# gepol_H3+_from-file.cpp test
add_Catch_test(
  NAME
    gepol_H3+_from-file
  LABELS
    gepol
    gepol_H3+_from-file
  DEPENDS
    gepol_H3+
  )

# gepol_C2H4_D2h.cpp test
add_Catch_test(
  NAME
    gepol_C2H4_D2h
  LABELS
    gepol
    gepol_C2H4_D2h
  )

# gepol_H3+_C2v.cpp test
add_Catch_test(
  NAME
    gepol_H3+_C2v
  LABELS
    gepol
    gepol_H3+_C2v
  )

# gepol_C6H6.cpp test
add_Catch_test(
  NAME
    gepol_C6H6
  LABELS
    gepol
    gepol_C6H6
  )

# gepol_CO2_symmetry.cpp test
add_Catch_test(
  NAME
    gepol_CO2_symmetry
  LABELS
    gepol
    gepol_CO2_symmetry
  )

# gepol_NH3.cpp test
add_Catch_test(
  NAME
    gepol_NH3
  LABELS
    gepol
    gepol_NH3
  )

# gepol_point.cpp test
add_Catch_test(
  NAME
    gepol_point
  LABELS
    gepol
    gepol_point
  )

# gepol_H3+.cpp test
add_Catch_test(
  NAME
    gepol_H3+
  LABELS
    gepol
    gepol_H3+
  )

# gepol_NH3_from-file.cpp test
add_Catch_test(
  NAME
    gepol_NH3_from-file
  LABELS
    gepol
    gepol_NH3_from-file
  DEPENDS
    gepol_NH3
  )
