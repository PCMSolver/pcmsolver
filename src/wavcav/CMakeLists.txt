set(headers_list "")
# Get headers with .h extension
get_h_headers(headers_list)
# If you want to remove some headers specify them explictly here
list(REMOVE_ITEM headers_list "test.h")

set(sources_list "")
# Get sources with .c extension
get_c_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "main.c")

if(ENABLE_TESTS)
	add_library(wavcav STATIC ${sources_list})
	set(LIBS_TO_MERGE wavcav ${LIBS_TO_MERGE} PARENT_SCOPE)
	install(TARGETS wavcav ARCHIVE DESTINATION lib)
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
