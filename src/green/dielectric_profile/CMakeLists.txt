# List of headers
list(APPEND headers_list Anisotropic.hpp MembraneTanh.hpp Metal.hpp OneLayerErf.hpp OneLayerTanh.hpp ProfileTypes.hpp Sharp.hpp Uniform.hpp Yukawa.hpp)

# List of sources
list(APPEND sources_list )

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${PROJECT_SOURCE_DIR}/src/green/dielectric_profile)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/green/dielectric_profile)
endforeach()
