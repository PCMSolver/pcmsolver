# List of headers
list(APPEND _headers_list Anisotropic.hpp MembraneTanh.hpp Metal.hpp OneLayerErf.hpp OneLayerTanh.hpp ProfileTypes.hpp Sharp.hpp Uniform.hpp Yukawa.hpp)

# List of sources
list(APPEND _sources_list )

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(_to_count "${_sources_list}" "${_headers_list}")
write_to_cloc_list("${_to_count}")

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
# Sets install directory for all the headers in the list
foreach(_header ${_headers_list})
   install(FILES ${_header} DESTINATION include/green/dielectric_profile)
endforeach()
