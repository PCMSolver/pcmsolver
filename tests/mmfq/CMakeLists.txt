target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/mmfq-ohno.cpp
  )

add_Catch_test(
  NAME
    mmfq-ohno
  LABELS
    mmfq
    ohno
  REFERENCE_FILES
    FQ-2xH2O.npy
  )
