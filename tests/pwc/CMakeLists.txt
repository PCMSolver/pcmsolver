set(list_of_libraries solver green bi_operators cavity wemFiles wemSolvers wavcav utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
# pwc_point test not yet done (wavelet cavity cannot generate a cavity from a single sphere) 
#    add_boosttest(pwc_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwc_NH3  "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwc_C6H6 "${list_of_libraries}" "${external_libraries}")
