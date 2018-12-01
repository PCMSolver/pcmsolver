# Copy input files to ${PROJECT_BINARY_DIR}/tests/standalone (aka ${CMAKE_CURRENT_BINARY_DIR})
file(
  COPY
    standalone.pcm
    standalone_atoms.pcm
    standalone_explicit.pcm
    biphenyl_gauss_failure.pcm
    pNA_SPD_failure.pcm
  DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}
  )

add_test(
  NAME
    pcmsolver_standalone
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone
  PROPERTIES
  LABELS "pcmsolver;standalone"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-atoms
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_atoms.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-atoms
  PROPERTIES
  LABELS "pcmsolver;standalone;atoms"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-explicit
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_explicit.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-explicit
  PROPERTIES
  LABELS "pcmsolver;standalone;explicit"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

# Demonstrates catching non-SPD S matrix failure
add_test(
  NAME
  pcmsolver_pNA-SPD-failure
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=pNA_SPD_failure.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_pNA-SPD-failure
  PROPERTIES
  LABELS "pcmsolver;standalone;SPD-failure"
  WILL_FAIL TRUE
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

# Demonstrates failure caught by Gauss' theorem test
add_test(
  NAME
    pcmsolver_biphenyl-gauss-failure
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=biphenyl_gauss_failure.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_biphenyl-gauss-failure
  PROPERTIES
  LABELS "pcmsolver;standalone;gauss-failure"
  WILL_FAIL TRUE
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
