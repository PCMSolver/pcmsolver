# List of headers
list(APPEND headers_list InputManager.hpp Meddle.hpp)

# List of sources
list(APPEND sources_list Meddle.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(interface STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(interface cavity solver utils)
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES interface)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(interface custom_boost)
endif()
install(TARGETS interface ARCHIVE DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/interface)
endforeach()
