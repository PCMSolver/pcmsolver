set(Doxygen_FIND_QUIETLY TRUE)
include(FindDoxygen)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile
    )
    add_custom_target(
        doxygen 
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
else()
	message("Doxygen missing. You won't be able to create docs.")
endif()
