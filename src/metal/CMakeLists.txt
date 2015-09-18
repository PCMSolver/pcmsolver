# List of sources
list(APPEND sources_list metal_sphere.F90)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}")
write_to_cloc_list("${to_count}")

foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
