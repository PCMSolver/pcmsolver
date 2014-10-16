set(headers_list "")
# Get headers with .hpp extension
get_hpp_headers(headers_list)
# If you want to remove some headers specify them explictly here
list(REMOVE_ITEM headers_list "")
# Sort alphabetically
list(SORT headers_list)
set(to_includer "${headers_list}")
# Write headers to Includer.hpp
file(APPEND "${PROJECT_BINARY_DIR}/include/Includer.hpp" "// Construction of boundary integral operators\n")
list(REMOVE_ITEM to_includer "")
write_to_Includer("${to_includer}")

set(sources_list "")
# Get sources with .cpp extension
get_cpp_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "")

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(bi_operators STATIC ${sources_list})
set(libs_to_merge bi_operators ${libs_to_merge} PARENT_SCOPE)
add_dependencies(bi_operators metal libtaylor)
if(BUILD_CUSTOM_BOOST)
	add_dependencies(bi_operators custom_boost)
endif()	
install(TARGETS bi_operators ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include)
