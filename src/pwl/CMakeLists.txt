set(headers_list "")
# Get headers with .h extension
get_h_headers(headers_list)
# If you want to remove some headers specify them explictly here
list(REMOVE_ITEM headers_list "")

set(sources_list "")
# Get sources with .c extension
get_c_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "main_pwl.c")

# Build static library
add_library(pwl STATIC ${sources_list})
set(libs_to_merge pwl ${libs_to_merge} PARENT_SCOPE)
install(TARGETS pwl ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include)
