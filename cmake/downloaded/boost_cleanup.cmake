# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

# Clean-up
add_custom_command(
    OUTPUT ${CUSTOM_BOOST_LOCATION}/boost.cleanedup
    if(NOT WIN32)
        # observed trouble deleting the directory on Windows + MinGW + multiple builds
        COMMAND ${CMAKE_COMMAND} -E rm -r ${BOOST_BUILD_DIR}
    endif()
    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_BOOST_LOCATION}/boost.cleanedup
    WORKING_DIRECTORY ${CUSTOM_BOOST_LOCATION}
    DEPENDS ${CUSTOM_BOOST_LOCATION}/boost.installed
    COMMENT "Clean-up Boost")
