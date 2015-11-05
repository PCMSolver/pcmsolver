set_property(GLOBAL PROPERTY TestSources)
set(TESTS_ROOT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(surface_function)
add_subdirectory(input)
add_subdirectory(numerical_quadrature)
add_subdirectory(gepol)
add_subdirectory(dielectric_profile)
add_subdirectory(green)
add_subdirectory(bi_operators)
add_subdirectory(cpcm)
add_subdirectory(iefpcm)

get_property(test_sources GLOBAL PROPERTY TestSources)

add_executable(unit_tests.x unit_tests.cpp ${test_sources})
target_link_libraries(unit_tests.x pcm getkw ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})

add_executable(C_host.x C_host.c)
target_link_libraries(C_host.x pcm getkw ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
set_target_properties(C_host.x PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests
                                          INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                          LINKER_LANGUAGE C)
add_test(NAME C_host COMMAND C_host.x)

if(ENABLE_FORTRAN_API)
    add_executable(Fortran_host.x Fortran_host.F90 Fortran_host-modules.F90)
    target_link_libraries(Fortran_host.x pcm getkw ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
    set_target_properties(Fortran_host.x PROPERTIES LINKER_LANGUAGE Fortran)
    add_test(NAME Fortran_host COMMAND Fortran_host.x)
endif()
