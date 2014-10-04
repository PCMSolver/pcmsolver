set(list_of_libraries solver green bi_operators cavity pedra utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES})
add_boosttest(iefpcm_gepol-point "${list_of_libraries}" "${external_libraries}")
add_boosttest(iefpcm_gepol-NH3 "${list_of_libraries}" "${external_libraries}")
add_boosttest(iefpcm_gepol-point_from-file "${list_of_libraries}" "${external_libraries}")  
add_boosttest(iefpcm_gepol-NH3_from-file "${list_of_libraries}" "${external_libraries}")  

if(DEVELOPMENT_CODE)
    set(list_of_libraries solver green bi_operators cavity tsless pedra utils)
    add_boosttest(iefpcm_tsless-point "${list_of_libraries}" "${external_libraries}")
    add_boosttest(iefpcm_tsless-NH3 "${list_of_libraries}" "${external_libraries}")
endif()

list(APPEND external_libraries ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
add_boosttest(iefpcm_symmetry "${list_of_libraries}" "${external_libraries}")  
