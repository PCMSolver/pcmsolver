set(list_of_libraries solver green bi_operators cavity wemFiles wemSolvers wavcav utils)
# wem_lin_point test not yet done 
#    add_boosttest(wem_lin_point "${list_of_libraries}" "${LIBS}")
add_boosttest(wem_lin_NH3 "${list_of_libraries}" "${LIBS}")
