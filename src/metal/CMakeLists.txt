# List of sources
list(APPEND sources_list metal_sphere.F90)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(metal STATIC ${sources_list})
# Specify dependencies for the library (if any)
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES metal)
install(TARGETS metal ARCHIVE DESTINATION lib)

