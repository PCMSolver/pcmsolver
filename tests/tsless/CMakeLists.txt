set(list_of_libraries cavity tsless pedra utils)
add_boosttest(tsless_point "${list_of_libraries}" "${LIBS}")
add_boosttest(tsless_NH3 "${list_of_libraries}" "${LIBS}")
