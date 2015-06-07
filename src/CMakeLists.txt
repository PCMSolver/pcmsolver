# List here all directories that contain *.h or *.hpp
# header files. This is needed in order to set the list
# of directories where CMake will go look for headers
set(header_dir_list
	${CMAKE_CURRENT_SOURCE_DIR}/cavity
	${CMAKE_CURRENT_SOURCE_DIR}/green
	${CMAKE_CURRENT_SOURCE_DIR}/green/dielectric_profile
	${CMAKE_CURRENT_SOURCE_DIR}/interface
    ${CMAKE_CURRENT_SOURCE_DIR}/bi_operators # Stands for Boundary Integral operators
	${CMAKE_CURRENT_SOURCE_DIR}/solver
	${CMAKE_CURRENT_SOURCE_DIR}/utils
    )
include_directories(${header_dir_list})

init_Includer()

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(bi_operators)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)
# This directory contains source file for executables
if(BUILD_STANDALONE)
   add_subdirectory(bin)
endif()

# We need to export header_dir_list one level up
set(header_dir_list ${header_dir_list} PARENT_SCOPE)

set(libs_to_merge ${libs_to_merge} PARENT_SCOPE)
