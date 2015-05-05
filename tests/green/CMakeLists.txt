set(list_of_libraries green bi_operators)

list(APPEND external_libraries ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
add_boosttest(green_vacuum "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_uniform_dielectric "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_ionic_liquid "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_anisotropic_liquid "${list_of_libraries}" "${external_libraries}")
add_boosttest(green_tanh_spherical_diffuse "${list_of_libraries}" "${external_libraries}")
