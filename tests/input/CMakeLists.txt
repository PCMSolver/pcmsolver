add_library(input-tests OBJECT
  input_cpcm.cpp
  input_diffuse.cpp
  input_multipoles.cpp
  input_restart.cpp
  )
set_target_properties(input-tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(input-tests custom_boost)
endif()

file(
  COPY
    diffuse.inp
    restart.inp
    cpcm.inp
    multipoles.inp
  DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}
  )

add_test(
  NAME
    parser_script-diffuse
  COMMAND
    ${PYTHON_EXECUTABLE} "-c" "from pcmsolver import parse_pcm_input; \
     parsed = parse_pcm_input('diffuse.inp', write_out=True)"
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(parser_script-diffuse
  PROPERTIES
  LABELS "parser;parser_script"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
add_Catch_test(input_diffuse "input;input_diffuse")
set_tests_properties(input_diffuse PROPERTIES DEPENDS parser_script-diffuse)

add_test(
  NAME
    parser_script-multipoles
  COMMAND
    ${PYTHON_EXECUTABLE} "-c" "from pcmsolver import parse_pcm_input; \
     parsed = parse_pcm_input('multipoles.inp', write_out=True)"
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(parser_script-multipoles
  PROPERTIES
  LABELS "parser;parser_script"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
add_Catch_test(input_multipoles "input;input_multipoles")
set_tests_properties(input_multipoles PROPERTIES DEPENDS parser_script-multipoles)

add_test(
  NAME
    parser_script-restart
  COMMAND
    ${PYTHON_EXECUTABLE} "-c" "from pcmsolver import parse_pcm_input; \
     parsed = parse_pcm_input('restart.inp', write_out=True)"
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(parser_script-restart
  PROPERTIES
  LABELS "parser;parser_script"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
add_Catch_test(input_restart "input;input_restart")
set_tests_properties(input_restart PROPERTIES DEPENDS parser_script-restart)

add_test(
  NAME
    parser_script-cpcm
  COMMAND
    ${PYTHON_EXECUTABLE} "-c" "from pcmsolver import parse_pcm_input; \
     parsed = parse_pcm_input('cpcm.inp', write_out=True)"
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(parser_script-cpcm
  PROPERTIES
  LABELS "parser;parser_script"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
add_Catch_test(input_cpcm "input;input_cpcm")
set_tests_properties(input_cpcm PROPERTIES DEPENDS parser_script-cpcm)
