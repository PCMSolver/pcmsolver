# List of sources
list(APPEND sources_list
tsless_cavity.f90
tsless_lapack.f90
tsless_precision.f90
tsless_symmetry.f90
tsless_weight_function.f90
)

add_library(tsless OBJECT ${sources_list})
set_target_properties(tsless
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
target_compile_definitions(tsless
  PUBLIC
    PCMSolver_EXPORTS
  )
