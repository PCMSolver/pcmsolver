# List of headers
list(APPEND headers_list Cubature.hpp GaussLegendre.hpp GaussSquare.hpp Quadrature.hpp WEMPCG.hpp WEMPGMRES.hpp WEMRHS.hpp)

# List of sources
list(APPEND sources_list GaussLegendre.cpp GaussSquare.cpp WEMPCG.cpp WEMPGMRES.cpp WEMRHS.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(wemSolvers STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(wemSolvers )
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES wemSolvers)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(wemSolvers custom_boost)
endif()
install(TARGETS wemSolvers ARCHIVE DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${PROJECT_SOURCE_DIR}/wemSolvers)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/wemSolvers)
endforeach()
