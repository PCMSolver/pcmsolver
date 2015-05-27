set(list_of_libraries)

list(APPEND external_libraries "")
add_boosttest(tanh_one_layer "${list_of_libraries}" "${external_libraries}")
