add_library(input-tests OBJECT input.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(input-tests custom_boost)
endif()

# Copy input files to ${PROJECT_BINARY_DIR}/tests/input (aka ${CMAKE_CURRENT_BINARY_DIR})
file(COPY diffuse.inp restart.inp tsless.inp wavelet.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# Tests pcmsolver.py
add_test(NAME parser_script-diffuse COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py diffuse.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-diffuse PROPERTIES LABELS "parser_script" DEPENDS generate-pcmsolver-py)
add_test(NAME parser_script-restart COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py restart.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-restart PROPERTIES LABELS "parser_script" DEPENDS generate-pcmsolver-py)
add_test(NAME parser_script-tsless  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py tsless.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-tsless  PROPERTIES LABELS "parser_script" DEPENDS generate-pcmsolver-py)
add_test(NAME parser_script-wavelet COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py wavelet.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-wavelet PROPERTIES LABELS "parser_script" DEPENDS generate-pcmsolver-py)

# input.cpp test
add_Catch_test(input "input")
set_tests_properties(input PROPERTIES DEPENDS parser_script-diffuse
                                      DEPENDS parser_script-restart
                                      DEPENDS parser_script-tsless
                                      DEPENDS parser_script-wavelet)
