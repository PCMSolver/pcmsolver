set(list_of_libraries "")
add_boosttest(numerical_quadrature_boost "${list_of_libraries}" "${LIBS}")
list(APPEND list_of_libraries cavity pedra utils)
add_boosttest(numerical_quadrature "${list_of_libraries}" "${LIBS}")
