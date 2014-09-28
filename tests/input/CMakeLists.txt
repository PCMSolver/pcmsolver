set(list_of_libraries utils getkw ${PYTHON_LIBRARIES})
add_boosttest(input "${list_of_libraries}" "${LIBS}")
add_reference(parsingTest.inp ${CMAKE_CURRENT_BINARY_DIR})
