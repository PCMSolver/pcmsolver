set(list_of_libraries solver green cavity pedra utils)
add_googletest(iefpcm_gepol-point "${list_of_libraries}" "${LIBS}")
add_googletest(iefpcm_gepol-NH3 "${list_of_libraries}" "${LIBS}")

if(DEVELOPMENT_CODE)
    set(list_of_libraries solver green cavity tsless pedra utils)
    add_googletest(iefpcm_tsless-point "${list_of_libraries}" "${LIBS}")
    add_googletest(iefpcm_tsless-NH3 "${list_of_libraries}" "${LIBS}")
endif()
