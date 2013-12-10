set(list_of_libraries solver green cavity pwl wem wavcav utils)
# The pwc_point test is not compiled as it does not work: Issue #547
#    add_googletest(pwc_point "${list_of_libraries}" "${LIBS}")
add_googletest(pwc_NH3 "${list_of_libraries}" "${LIBS}")
