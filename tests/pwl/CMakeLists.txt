set(list_of_libraries solver green bi_operators cavity pwl wem wavcav utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
# The pwl_point test is not compiled as it does not work: Issue #547
#    add_boosttest(pwl_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(pwl_NH3 "${list_of_libraries}" "${external_libraries}")
