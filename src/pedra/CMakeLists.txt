# List of sources
list(APPEND sources_list
       pedra_cavity.F90
       pedra_cavity_interface.F90
       pedra_dblas.f90
       pedra_dlapack.f90
       pedra_precision.f90
       pedra_print.f90
       pedra_symmetry.f90
       pedra_utils.f90
    )

add_library(pedra OBJECT ${sources_list})
set_target_properties(pedra PROPERTIES POSITION_INDEPENDENT_CODE 1 )
