add_library(input-tests OBJECT input_diffuse.cpp input_restart.cpp input_tsless.cpp input_wavelet.cpp input_multipoles.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(input-tests custom_boost)
endif()

# Copy input files to ${PROJECT_BINARY_DIR}/tests/input (aka ${CMAKE_CURRENT_BINARY_DIR})
file(COPY diffuse.inp restart.inp
          tsless.inp wavelet.inp
          multipoles.inp
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME parser_script-diffuse COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py diffuse.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-diffuse PROPERTIES LABELS "parser;parser_script" DEPENDS generate-pcmsolver-py)
add_Catch_test(input_diffuse "input;input_diffuse")
set_tests_properties(input_diffuse PROPERTIES DEPENDS parser_script-diffuse)

add_test(NAME parser_script-multipoles COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py multipoles.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-multipoles PROPERTIES LABELS "parser;parser_script" DEPENDS generate-pcmsolver-py)
add_Catch_test(input_multipoles "input;input_multipoles")
set_tests_properties(input_multipoles PROPERTIES DEPENDS parser_script-multipoles)

add_test(NAME parser_script-restart COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py restart.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-restart PROPERTIES LABELS "parser;parser_script" DEPENDS generate-pcmsolver-py)
add_Catch_test(input_restart "input;input_restart")
set_tests_properties(input_restart PROPERTIES DEPENDS parser_script-restart)

add_test(NAME parser_script-tsless  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py tsless.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-tsless  PROPERTIES LABELS "parser;parser_script" DEPENDS generate-pcmsolver-py)
add_Catch_test(input_tsless "input;input_tsless")
set_tests_properties(input_tsless PROPERTIES DEPENDS parser_script-tsless)

add_test(NAME parser_script-wavelet COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py wavelet.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parser_script-wavelet PROPERTIES LABELS "parser;parser_script" DEPENDS generate-pcmsolver-py)
add_Catch_test(input_wavelet "input;input_wavelet")
set_tests_properties(input_wavelet PROPERTIES DEPENDS parser_script-wavelet)
