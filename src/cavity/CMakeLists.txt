# List of headers
list(APPEND headers_list Cavity.hpp Element.hpp GePolCavity.hpp RegisterCavityToFactory.hpp RestartCavity.hpp)

# List of sources
list(APPEND sources_list Cavity.cpp Element.cpp GePolCavity.cpp RestartCavity.cpp)

# Update bar chart
if(BUILD_CHARTS)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR} CXX)
endif()

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/cavity)
endforeach()
