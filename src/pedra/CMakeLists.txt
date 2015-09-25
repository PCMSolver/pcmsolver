# List of sources
list(APPEND sources_list pedra_cavity.F90 pedra_cavity_interface.F90 pedra_dblas.F90 pedra_dlapack.F90 pedra_precision.F90 pedra_print.F90 pedra_symmetry.F90 pedra_utils.F90)

# Update bar chart
if(BUILD_DOCUMENTATION)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR} Fortran)
endif()

foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()

