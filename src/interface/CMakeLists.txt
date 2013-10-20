set(interface_SOURCES 
        Interface.cpp
)

if(ENABLE_TESTS)
	add_library(interface ${interface_SOURCES})
	add_dependencies(interface utils solver)
	install(TARGETS interface ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${interface_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(FILES Interface.hpp DESTINATION include)
