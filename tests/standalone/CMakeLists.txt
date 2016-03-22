# Add test using the standalone executable
add_test(NAME pcmsolver_standalone COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py -x standalone.pcm)
set_tests_properties(pcmsolver_standalone PROPERTIES LABELS "pcmsolver;standalone")
