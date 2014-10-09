set(headers_list "")
# Get headers with .hpp extension
get_hpp_headers(headers_list)
# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
	list(REMOVE_ITEM headers_list "")
else()
	list(REMOVE_ITEM headers_list "WEMSolver.hpp")
endif()	
# Sort alphabetically
list(SORT headers_list)
# Write headers to Includer.hpp
file(APPEND "${PROJECT_BINARY_DIR}/include/Includer.hpp" "// Solvers\n")
write_to_Includer("${headers_list}" "")

set(sources_list "")
# Get sources with .cpp extension
get_cpp_sources(sources_list)
# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
	list(REMOVE_ITEM sources_list "")
else()
	list(REMOVE_ITEM sources_list "WEMSolver.cpp")
endif()	

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(solver STATIC ${sources_list})
set(libs_to_merge solver ${libs_to_merge} PARENT_SCOPE)
if(DEVELOPMENT_CODE)
	add_dependencies(solver cavity green wem pwl)
else()
	add_dependencies(solver cavity green)
endif()
install(TARGETS solver ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include)
