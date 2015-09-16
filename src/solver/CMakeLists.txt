# List of headers
list(APPEND headers_list CPCMSolver.hpp IEFSolver.hpp PCMSolver.hpp RegisterSolverToFactory.hpp)

# List of sources
list(APPEND sources_list CPCMSolver.cpp IEFSolver.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(solver STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(solver cavity)
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES solver)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(solver custom_boost)
endif()
install(TARGETS solver ARCHIVE DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/solver)
endforeach()
