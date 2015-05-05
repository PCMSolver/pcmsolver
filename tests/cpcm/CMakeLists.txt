set(list_of_libraries solver green bi_operators cavity pedra utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
add_boosttest(cpcm_gepol-point "${list_of_libraries}" "${external_libraries}")
add_boosttest(cpcm_gepol-NH3 "${list_of_libraries}" "${external_libraries}")
add_boosttest(cpcm_gepol-point_from-file "${list_of_libraries}" "${external_libraries}")
add_boosttest(cpcm_gepol-NH3_from-file "${list_of_libraries}" "${external_libraries}")
add_boosttest(cpcm_symmetry "${list_of_libraries}" "${external_libraries}")

