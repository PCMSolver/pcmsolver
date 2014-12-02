# Executables for testing separate from Boost unit test framework or for profiling purposes
# The executables will be created in ${PROJECT_BINARY_DIR}/src/bin

#add_executable(debug_wavcav.x debug_wavcav.cpp)
#target_link_libraries(debug_wavcav.x solver green cavity wavcav pwl wem utils)

add_executable(profile_pwc.x profile_pwc.cpp)
target_link_libraries(profile_pwc.x solver green bi_operators cavity wemFiles wemSolvers wavcav utils ${ZLIB_LIBRARY})

add_executable(profile_pwl.x profile_pwl.cpp)
target_link_libraries(profile_pwl.x solver green bi_operators cavity wemFiles wemSolvers wavcav utils ${ZLIB_LIBRARY})
