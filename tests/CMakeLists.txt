set_property(GLOBAL PROPERTY TestSources)
# Math utilities and helper classes
add_subdirectory(surface_function)
#add_subdirectory(input)
#add_subdirectory(numerical_quadrature)
# Cavity generators
#add_subdirectory(gepol)
# Green's functions
#add_subdirectory(dielectric_profile)
#add_subdirectory(green)
#add_subdirectory(bi_operators)
# Solvers
#add_subdirectory(cpcm)
#add_subdirectory(iefpcm)

get_property(test_sources GLOBAL PROPERTY TestSources)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
add_executable(unit_tests.x unit_tests.cpp ${test_sources})
target_link_libraries(unit_tests.x pcm)
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

