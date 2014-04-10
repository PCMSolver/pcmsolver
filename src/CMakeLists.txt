set(header_dir_list 
	${CMAKE_CURRENT_SOURCE_DIR}/cavity 
	${CMAKE_CURRENT_SOURCE_DIR}/green
	${CMAKE_CURRENT_SOURCE_DIR}/interface
	${CMAKE_CURRENT_SOURCE_DIR}/pwl
	${CMAKE_CURRENT_SOURCE_DIR}/solver
	${CMAKE_CURRENT_SOURCE_DIR}/utils
	${CMAKE_CURRENT_SOURCE_DIR}/wavcav
	${CMAKE_CURRENT_SOURCE_DIR}/wem
   )
include_directories(${header_dir_list})

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(pedra)
add_subdirectory(pwl)
add_subdirectory(solver)
add_subdirectory(utils)
add_subdirectory(wem)
if(DEVELOPMENT_CODE)
	add_subdirectory(tsless)
	add_subdirectory(wavcav)
endif()

# We need to export header_dir_list one level up
set(header_dir_list ${header_dir_list} PARENT_SCOPE)

set(SOURCES ${SOURCES} PARENT_SCOPE)
