# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

# Clean-up
if(WIN32)
  # observed trouble deleting the directory on Windows + MinGW + multiple builds
  add_custom_command(
    OUTPUT ${CUSTOM_BOOST_LOCATION}/boost.cleanedup
    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_BOOST_LOCATION}/boost.cleanedup
    WORKING_DIRECTORY ${CUSTOM_BOOST_LOCATION}
    DEPENDS ${CUSTOM_BOOST_LOCATION}/boost.installed
    COMMENT "Clean-up Boost")
else()
  add_custom_command(
    OUTPUT ${CUSTOM_BOOST_LOCATION}/boost.cleanedup
    COMMAND ${CMAKE_COMMAND} -E rm -r ${BOOST_BUILD_DIR}
    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_BOOST_LOCATION}/boost.cleanedup
    WORKING_DIRECTORY ${CUSTOM_BOOST_LOCATION}
    DEPENDS ${CUSTOM_BOOST_LOCATION}/boost.installed
    COMMENT "Clean-up Boost")
endif()
