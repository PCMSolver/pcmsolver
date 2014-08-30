set(Doxygen_FIND_QUIETLY TRUE)
include(FindDoxygen)

macro(update_bar_charts)
	execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/counter.py")
	file(GLOB scripts "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/*.py")
        foreach(script ${scripts})                                                               
        	# Update bar charts
                execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${script}" 
                                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx)
        endforeach()
endmacro()

if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile)
    message(STATUS "Updating bar charts...")
    update_bar_charts()
    message(STATUS "    ...Done!")
    add_custom_target(
        doc 
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    add_custom_target(
	update_gh-pages
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/update_gh-pages.py"
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
else()
	message(STATUS "Doxygen missing. You won't be able to create docs.")
endif()
