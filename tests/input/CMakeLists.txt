set(list_of_libraries utils getkw ${PYTHON_LIBRARIES})
add_boosttest(input "${list_of_libraries}" "${LIBS}")
add_reference(gepol.inp ${CMAKE_CURRENT_BINARY_DIR})
add_reference(restart.inp ${CMAKE_CURRENT_BINARY_DIR})
