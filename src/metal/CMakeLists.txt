# List of sources
list(APPEND sources_list metal_sphere.F90)

add_library(metal OBJECT ${sources_list})
set_target_properties(metal PROPERTIES POSITION_INDEPENDENT_CODE 1 )
