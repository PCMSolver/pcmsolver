# Math utilities and helper classes
add_subdirectory(block_diagonal)
add_subdirectory(surface_function)
add_subdirectory(input)
# Cavity generators
add_subdirectory(gepol)
if(DEVELOPMENT_CODE)
	add_subdirectory(tsless)
	add_subdirectory(wavelet_cavity)
endif()
# Green's functions
add_subdirectory(green)
# Solvers
add_subdirectory(cpcm)
add_subdirectory(iefpcm)
if(DEVELOPMENT_CODE)
	#	add_subdirectory(pwc)
	#	add_subdirectory(pwl)
	add_subdirectory(wem_lin)
endif()

