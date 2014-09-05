set(list_of_libraries solver green operators_diagonal cavity pwl wem wavcav utils)
# The pwc_point test is not compiled as it does not work: Issue #547
#    add_boosttest(pwc_point "${list_of_libraries}" "${LIBS}")
add_boosttest(pwc_NH3 "${list_of_libraries}" "${LIBS}")
