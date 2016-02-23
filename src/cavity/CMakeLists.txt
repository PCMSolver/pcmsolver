# List of headers
list(APPEND headers_list Cavity.hpp CavityData.hpp Element.hpp GePolCavity.hpp RegisterCavityToFactory.hpp RestartCavity.hpp)

# List of sources
list(APPEND sources_list Cavity.cpp Element.cpp GePolCavity.cpp RestartCavity.cpp)

add_library(cavity OBJECT ${sources_list} ${headers_list})
set_target_properties(cavity PROPERTIES POSITION_INDEPENDENT_CODE 1 INTERPROCEDURAL_OPTIMIZATION 1)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/cavity)
endforeach()
