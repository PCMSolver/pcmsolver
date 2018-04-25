target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/input_cpcm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input_diffuse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input_multipoles.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input_restart.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input_sharp.cpp
  )

file(
  COPY
    cpcm.inp
    diffuse.inp
    multipoles.inp
    restart.inp
    sharp.inp
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
add_Catch_test(
  NAME
    input_diffuse
  LABELS
    input
    input_diffuse
  DEPENDS
    parser_script-diffuse
  )

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
add_Catch_test(
  NAME
    input_multipoles
  LABELS
    input
    input_multipoles
  DEPENDS
    parser_script-multipoles
  )

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
add_Catch_test(
  NAME
    input_restart
  LABELS
    input
    input_restart
  DEPENDS
    parser_script-restart
  )


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
add_Catch_test(
  NAME
    input_cpcm
  LABELS
    input
    input_cpcm
  DEPENDS
    parser_script-cpcm
  )

add_test(
  NAME
    parser_script-sharp
  COMMAND
    ${PYTHON_EXECUTABLE} "-c" "from pcmsolver import parse_pcm_input; \
     parsed = parse_pcm_input('sharp.inp', write_out=True)"
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(parser_script-sharp
  PROPERTIES
  LABELS "parser;parser_script"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
add_Catch_test(
  NAME
    input_sharp
  LABELS
    input
    input_sharp
  DEPENDS
    parser_script-sharp
  )
