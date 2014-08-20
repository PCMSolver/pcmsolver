set(Doxygen_FIND_QUIETLY TRUE)
include(FindDoxygen)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile
    )
    add_custom_target(
        doc 
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    add_custom_target(
	update_gh-pages
	COMMAND python ${PROJECT_BINARY_DIR}/bin/update_gh-pages.py
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	)
else()
	message(STATUS "Doxygen missing. You won't be able to create docs.")
endif()
