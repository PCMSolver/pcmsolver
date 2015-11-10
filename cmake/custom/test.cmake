option(ENABLE_TESTS "Enable PCMSolver unit tests" ON)

macro(add_Catch_test _name _labels)
  # _labels is not a list, it's a string... Transform it into a list
  set(labels)
  string(REPLACE ";" " " _labels "${_labels}")
  foreach(_label "${_labels}")
    list(APPEND labels ${_label})
  endforeach()
  unset(_labels)

  # This is the unit tests runner
  set(RUNNER ${PROJECT_BINARY_DIR}/tests/unit_tests${EXE})

  add_test(NAME ${_name} COMMAND ${RUNNER} [${_name}])

  if(labels)
    set_tests_properties(${_name} PROPERTIES LABELS "${labels}")
  endif()
endmacro()

if(ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(tests) # This must come last!!
endif()
