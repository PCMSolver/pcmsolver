set(list_of_libraries solver green bi_operators cavity wemFiles wemSolvers wavcav utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
# The pwl_point test is not compiled as it does not work: Issue #547
#    add_boosttest(pwl_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwl_NH3 "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwl_C6H6 "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwl_secondKind_C6H6 "${list_of_libraries}" "${external_libraries}")
