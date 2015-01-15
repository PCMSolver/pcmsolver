# Executables for testing separate from Boost unit test framework or for profiling purposes
# The executables will be created in ${PROJECT_BINARY_DIR}/src/bin

#add_executable(debug_wavcav.x debug_wavcav.cpp)
#target_link_libraries(debug_wavcav.x solver green cavity wavcav pwl wem utils)

add_executable(profile_pwc.x profile_pwc.cpp)
target_link_libraries(profile_pwc.x solver green bi_operators cavity wemFiles wemSolvers wavcav utils ${ZLIB_LIBRARY})

add_executable(profile_pwl.x profile_pwl.cpp)
target_link_libraries(profile_pwl.x solver green bi_operators cavity wemFiles wemSolvers wavcav utils ${ZLIB_LIBRARY})

include(TestingMacros) # imports add_reference

add_executable(benzene_wavelet.x benzene_wavelet.cpp)
target_link_libraries(benzene_wavelet.x solver green bi_operators cavity wemFiles wemSolvers wavcav utils ${ZLIB_LIBRARY})
add_reference(${CMAKE_SOURCE_DIR}/examples/benzene2.dat ${CMAKE_CURRENT_BINARY_DIR})
add_reference(${CMAKE_SOURCE_DIR}/examples/sphere2.dat ${CMAKE_CURRENT_BINARY_DIR})
