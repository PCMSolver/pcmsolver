# List of sources
list(APPEND sources_list tsless_cavity.F90 tsless_lapack.F90 tsless_precision.F90 tsless_symmetry.F90 tsless_weight_function.F90)

add_library(tsless OBJECT ${sources_list})
set_target_properties(tsless PROPERTIES POSITION_INDEPENDENT_CODE 1 )
