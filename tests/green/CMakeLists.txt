set(list_of_libraries green)
add_googletest(vacuum "${list_of_libraries}" "${LIBS}")
add_googletest(uniform_dielectric "${list_of_libraries}" "${LIBS}")
