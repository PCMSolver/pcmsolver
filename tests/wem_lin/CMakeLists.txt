set(list_of_libraries solver green bi_operators cavity wemFiles wemSolvers wavcav utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
# wem_lin_point test not yet done 
#    add_boosttest(wem_lin_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(wem_lin_NH3 "${list_of_libraries}" "${external_libraries}")
add_boosttest(wem_con_NH3 "${list_of_libraries}" "${external_libraries}")
