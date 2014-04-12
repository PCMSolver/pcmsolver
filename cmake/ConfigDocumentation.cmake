set(Doxygen_FIND_QUIETLY TRUE)
include(FindDoxygen)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile
    )
    configure_file(
      ${CMAKE_SOURCE_DIR}/doc/eigendoxy_header.html.in
      ${PROJECT_BINARY_DIR}/doc/eigendoxy_header.html
    )
    
    configure_file(
      ${CMAKE_SOURCE_DIR}/doc/eigendoxy_footer.html.in
      ${PROJECT_BINARY_DIR}/doc/eigendoxy_footer.html
    )
    
    configure_file(
      ${CMAKE_SOURCE_DIR}/doc/eigendoxy_layout.xml.in
      ${PROJECT_BINARY_DIR}/doc/eigendoxy_layout.xml
    )
    file(COPY ${CMAKE_SOURCE_DIR}/doc/eigen_navtree_hacks.js 
      DESTINATION ${PROJECT_BINARY_DIR}/doc/_doxygen/html
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE 
                       GROUP_READ 
                       WORLD_READ
        )
    add_custom_target(
        doxygen
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
else()
	message("Doxygen missing. You won't be able to create docs.")
endif()
