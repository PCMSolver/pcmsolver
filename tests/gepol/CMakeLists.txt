set(list_of_libraries cavity pedra utils)
# Order is important here, as the former two tests will produce the .npz files 
# to be read by the latter two tests.
add_boosttest(gepol_point "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_NH3 "${list_of_libraries}" "${LIBS}") 
add_boosttest(gepol_point_from-file "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_NH3_from-file "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_H3+ "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_H3+_from-file "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_point_symmetry "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_CO2_symmetry "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_CH3+_Cs "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_H3+_C2v "${list_of_libraries}" "${LIBS}")
add_boosttest(gepol_C2H4_D2h "${list_of_libraries}" "${LIBS}")
