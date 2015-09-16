# List of headers
list(APPEND headers_list Cavity.hpp Element.hpp GePolCavity.hpp RestartCavity.hpp)

# List of sources
list(APPEND sources_list Cavity.cpp Element.cpp GePolCavity.cpp RestartCavity.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(cavity STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cavity pedra)
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES cavity)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cavity custom_boost)
endif()
install(TARGETS cavity ARCHIVE DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${PROJECT_SOURCE_DIR}/src/cavity)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/cavity)
endforeach()
