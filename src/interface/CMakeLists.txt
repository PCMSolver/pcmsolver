set(headers_list "")
# Get headers with .hpp extension
get_hpp_headers(headers_list)
# If you want to remove some headers specify them explictly here
list(REMOVE_ITEM headers_list "")
# Sort alphabetically
#list(SORT headers_list)
# Write headers to Includer.hpp
#file(APPEND "${PROJECT_BINARY_DIR}/include/Includer.hpp" "// Interface\n")
#write_to_Includer("InputManager.hpp")

set(sources_list "")
# Get sources with .cpp extension
get_cpp_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "")

if(ENABLE_TESTS)
	add_library(interface ${sources_list})
	add_dependencies(interface utils solver)
	install(TARGETS interface ARCHIVE DESTINATION lib)
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
