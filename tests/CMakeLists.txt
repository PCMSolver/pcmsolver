# Math utilities and helper classes
add_subdirectory(block_diagonal)
add_subdirectory(surface_function)
add_subdirectory(input)
add_subdirectory(numerical_quadrature)
# Cavity generators
add_subdirectory(gepol)
if(DEVELOPMENT_CODE)
	add_subdirectory(tsless)
	add_subdirectory(wavelet_cavity)
endif()
# Green's functions
add_subdirectory(green)
add_subdirectory(bi_operators)
# Solvers
add_subdirectory(cpcm)
add_subdirectory(iefpcm)
if(DEVELOPMENT_CODE)
	add_subdirectory(pwc)
	add_subdirectory(pwl)
endif()

