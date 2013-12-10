set(list_of_libraries cavity pedra utils)
add_googletest(gepol_point "${list_of_libraries}" "${LIBS}")
add_googletest(gepol_NH3 "${list_of_libraries}" "${LIBS}") 
