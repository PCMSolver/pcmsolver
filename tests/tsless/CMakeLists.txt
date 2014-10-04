set(list_of_libraries cavity tsless pedra utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
add_boosttest(tsless_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(tsless_NH3 "${list_of_libraries}" "${external_libraries}")
