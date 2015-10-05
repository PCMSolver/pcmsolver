# List of sources
list(APPEND sources_list tsless_cavity.F90 tsless_lapack.F90 tsless_precision.F90 tsless_symmetry.F90 tsless_weight_function.F90)

# Update bar chart
if(BUILD_CHARTS)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR} Fortran)
endif()

foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
