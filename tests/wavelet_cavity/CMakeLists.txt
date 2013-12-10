if(DEVELOPMENT_CODE)
    set(list_of_libraries solver green cavity wem pwl wavcav utils)
    # The wavelet_cavity_point test is not compiled as it does not work: Issue #547
    #    add_googletest(wavelet_cavity_point "${list_of_libraries}" "${LIBS}")
    add_googletest(wavelet_cavity_NH3 "${list_of_libraries}" "${LIBS}")
endif()
