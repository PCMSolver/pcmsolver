set(list_of_libraries solver green bi_operators cavity pwl wem wavcav utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
# The pwc_point test is not compiled as it does not work: Issue #547
#    add_boosttest(pwc_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwc_NH3 "${list_of_libraries}" "${external_libraries}")
