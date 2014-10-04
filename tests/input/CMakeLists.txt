set(list_of_libraries utils getkw)

list(APPEND external_libraries ${PYTHON_LIBRARIES})
add_boosttest(input "${list_of_libraries}" "${external_libraries}")
add_reference(tsless.inp ${CMAKE_CURRENT_BINARY_DIR})
add_reference(restart.inp ${CMAKE_CURRENT_BINARY_DIR})
add_reference(wavelet.inp ${CMAKE_CURRENT_BINARY_DIR})
