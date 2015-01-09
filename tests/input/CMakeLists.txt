set(list_of_libraries utils getkw)

add_boosttest(input "${list_of_libraries}")
# Run pcmsolver.py on tsless.inp and copy @tsless.inp to the build directory
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/bin/pcmsolver.py tsless.inp
	        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_reference(@tsless.inp ${CMAKE_CURRENT_BINARY_DIR})
# Run pcmsolver.py on restart.inp and copy @restart.inp to the build directory
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/bin/pcmsolver.py restart.inp
	        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_reference(@restart.inp ${CMAKE_CURRENT_BINARY_DIR})
# Run pcmsolver.py on wavelet.inp and copy @wavelet.inp to the build directory
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/bin/pcmsolver.py wavelet.inp
	        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_reference(@wavelet.inp ${CMAKE_CURRENT_BINARY_DIR})
