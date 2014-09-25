set(headers_list "")
# Get headers with .h extension
get_h_headers(headers_list)
# If you want to remove some headers specify them explictly here
list(REMOVE_ITEM headers_list "")

set(sources_list "")
# Get sources with .c extension
get_c_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "data_bnz.c" "data.c" "data_pd.c" "init_exterior.c" "main2.c" "main.c" "matlab_c.c" "modify.c" "pot.c" "sphere.c")

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Override global compiler flags to avoid choking stdout with warnings
set_wem_compiler_flags()
# Build static library
add_library(wem STATIC ${sources_list})
set(libs_to_merge wem ${libs_to_merge} PARENT_SCOPE)
install(TARGETS wem ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include)
