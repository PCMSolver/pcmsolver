set(sources_list "")
# Get sources with Fortran extension
get_fortran_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "")

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(metal STATIC ${sources_list})
set(libs_to_merge metal ${libs_to_merge} PARENT_SCOPE)
install(TARGETS metal ARCHIVE DESTINATION lib)

