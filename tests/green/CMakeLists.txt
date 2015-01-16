set(list_of_libraries green bi_operators)

list(APPEND external_libraries)
add_boosttest(green_vacuum "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_uniform_dielectric "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_ionic_liquid "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_anisotropic_liquid "${list_of_libraries}" "${external_libraries}")
