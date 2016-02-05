# List of sources
list(APPEND sources_list pedra_cavity.F90 pedra_cavity_interface.F90 pedra_dblas.F90 pedra_dlapack.F90 pedra_precision.F90 pedra_print.F90 pedra_symmetry.F90 pedra_utils.F90)

add_library(pedra OBJECT ${sources_list})
set_target_properties(pedra PROPERTIES POSITION_INDEPENDENT_CODE 1 INTERPROCEDURAL_OPTIMIZATION 1)
