set(list_of_libraries solver green bi_operators cavity wemFiles wemSolvers wavcav utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
# The wavelet_cavity_point test is not compiled as it does not work: Issue #547
#    add_boosttest(wavelet_cavity_point "${list_of_libraries}" "${external_libraries}")
add_boosttest(wavelet_cavity_NH3 "${list_of_libraries}" "${external_libraries}")
add_boosttest(wavelet_cavity_C6H6 "${list_of_libraries}" "${external_libraries}")
