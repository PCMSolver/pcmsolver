# Copy input files to ${PROJECT_BINARY_DIR}/tests/standalone (aka ${CMAKE_CURRENT_BINARY_DIR})
file(COPY standalone.pcm standalone_atoms.pcm standalone_explicit.pcm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Add test using the standalone executable
add_test(NAME pcmsolver_standalone
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py -x standalone.pcm
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(pcmsolver_standalone PROPERTIES LABELS "pcmsolver;standalone"
                                                     DEPENDS generate-pcmsolver-py)

add_test(NAME pcmsolver_standalone-atoms
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py -x standalone_atoms.pcm
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(pcmsolver_standalone-atoms PROPERTIES LABELS "pcmsolver;standalone;atoms"
                                                     DEPENDS generate-pcmsolver-py)

add_test(NAME pcmsolver_standalone-explicit
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py -x standalone_explicit.pcm
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(pcmsolver_standalone-explicit PROPERTIES LABELS "pcmsolver;standalone;explicit"
                                                     DEPENDS generate-pcmsolver-py)
