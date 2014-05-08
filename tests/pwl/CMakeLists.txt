set(list_of_libraries solver green operators_diagonal cavity pwl wem wavcav utils)
# The pwl_point test is not compiled as it does not work: Issue #547
#    add_boosttest(pwl_point "${list_of_libraries}" "${LIBS}")
add_boosttest(pwl_NH3 "${list_of_libraries}" "${LIBS}")
