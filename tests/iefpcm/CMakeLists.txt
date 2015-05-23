set(list_of_libraries solver green cavity pedra utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
add_boosttest(iefpcm_gepol-point "${list_of_libraries}" "${external_libraries}")
#add_boosttest(iefpcm_diffuse-gepol-point "${list_of_libraries}" "${external_libraries}")
add_boosttest(iefpcm_gepol-NH3 "${list_of_libraries}" "${external_libraries}")
add_boosttest(iefpcm_gepol-point_from-file "${list_of_libraries}" "${external_libraries}")
add_boosttest(iefpcm_gepol-NH3_from-file "${list_of_libraries}" "${external_libraries}")
add_boosttest(iefpcm_symmetry "${list_of_libraries}" "${external_libraries}")
