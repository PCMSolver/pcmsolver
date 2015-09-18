# List of sources
list(APPEND sources_list pedra_cavity.F90 pedra_cavity_interface.F90 pedra_dblas.F90 pedra_dlapack.F90 pedra_precision.F90 pedra_print.F90 pedra_symmetry.F90 pedra_utils.F90)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}")
write_to_cloc_list("${to_count}")

foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()

