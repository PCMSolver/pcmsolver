set(list_of_libraries green)
add_googletest(green_vacuum "${list_of_libraries}" "${LIBS}")
add_googletest(green_uniform_dielectric "${list_of_libraries}" "${LIBS}")
add_googletest(green_ionic_liquid "${list_of_libraries}" "${LIBS}")
