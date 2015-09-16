# List of headers
list(APPEND headers_list CollocationIntegrator.hpp IntegratorForward.hpp IntegratorHelperFunctions.hpp IntegratorTypes.hpp NumericalIntegrator.hpp PurisimaIntegrator.hpp)

# List of sources
list(APPEND sources_list )

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${PROJECT_SOURCE_DIR}/src/bi_operators)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/bi_operators)
endforeach()
