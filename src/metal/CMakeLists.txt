# List of sources
list(APPEND sources_list metal_sphere.F90)

# Update bar chart
if(BUILD_DOCUMENTATION)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR})
endif()

foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
