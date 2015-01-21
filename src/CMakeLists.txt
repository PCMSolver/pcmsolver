# List here all directories that contain *.h or *.hpp
# header files. This is needed in order to set the list
# of directories where CMake will go look for headers
set(header_dir_list 
	${CMAKE_CURRENT_SOURCE_DIR}/cavity 
	${CMAKE_CURRENT_SOURCE_DIR}/green
	${CMAKE_CURRENT_SOURCE_DIR}/interface
	${CMAKE_CURRENT_SOURCE_DIR}/bi_operators # Stands for Boundary Integral operators
	${CMAKE_CURRENT_SOURCE_DIR}/logger
	${CMAKE_CURRENT_SOURCE_DIR}/solver
	${CMAKE_CURRENT_SOURCE_DIR}/utils)
if(DEVELOPMENT_CODE)
        set(header_dir_list_DEVELOPMENT            
        	${CMAKE_CURRENT_SOURCE_DIR}/pwl
        	${CMAKE_CURRENT_SOURCE_DIR}/wavcav
        	${CMAKE_CURRENT_SOURCE_DIR}/wem)
	list(APPEND header_dir_list "${header_dir_list_DEVELOPMENT}")
	list(SORT header_dir_list)
endif()	
include_directories(${header_dir_list})

init_Includer()

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(bi_operators)
add_subdirectory(logger)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)
if(DEVELOPMENT_CODE)
	# The alphabetic order rule is here violated as the
	# wem directory contains some files also needed by
	# wavcav and pwl...
	add_subdirectory(wem)
	add_subdirectory(pwl)
	add_subdirectory(tsless)
	add_subdirectory(wavcav)
endif()
# This directory contains source file for executables
add_subdirectory(bin)

# We need to export header_dir_list one level up
set(header_dir_list ${header_dir_list} PARENT_SCOPE)

set(libs_to_merge ${libs_to_merge} PARENT_SCOPE)
