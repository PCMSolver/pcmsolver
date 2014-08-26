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

if(ENABLE_TESTS)
	add_library(wem STATIC ${sources_list})
	set(LIBS_TO_MERGE wem ${LIBS_TO_MERGE} PARENT_SCOPE)
	install(TARGETS wem ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${sources_list})
	get_filename_component(sourcepath ${source} ABSOLUTE)
	list(APPEND local_sources ${sourcepath})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include)
