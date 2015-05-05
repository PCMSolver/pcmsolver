set(list_of_libraries "")
list(APPEND list_of_libraries cavity pedra utils)

list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
add_boosttest(numerical_quadrature "${list_of_libraries}" "${external_libraries}")
add_reference(molecule.npy ${CMAKE_CURRENT_BINARY_DIR})
add_reference(molecule_1r.npy ${CMAKE_CURRENT_BINARY_DIR})
