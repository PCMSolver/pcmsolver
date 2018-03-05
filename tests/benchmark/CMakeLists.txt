# Copy input files to ${PROJECT_BINARY_DIR}/tests/standalone (aka ${CMAKE_CURRENT_BINARY_DIR})
file(
  COPY
    standalone_bubble.pcm
    standalone_water-in-oil.pcm
    standalone_water-fake.pcm
  DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}
  )

add_test(
  NAME
    pcmsolver_standalone-bubble
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_bubble.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-bubble
  PROPERTIES
  LABELS "pcmsolver;standalone;bubble"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-water-in-oil
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_water-in-oil.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-water-in-oil
  PROPERTIES
  LABELS "pcmsolver;standalone;water-in-oil"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-water-fake
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_water-fake.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-water-fake
  PROPERTIES
  LABELS "pcmsolver;standalone;water-fake"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
