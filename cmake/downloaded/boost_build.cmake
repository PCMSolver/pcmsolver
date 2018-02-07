# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

# Build Boost
# This is not Windows-friendly!
add_custom_command(
    OUTPUT ${CUSTOM_BOOST_LOCATION}/boost.built
    COMMAND ./b2 toolset=${toolset} variant=${type} link=static cxxflags=-fPIC
    threading=multi --user-config=user-config.jam 1> ${CUSTOM_BOOST_LOCATION}/boost.built.log 2> ${CUSTOM_BOOST_LOCATION}/boost.built.err
    COMMAND ${CMAKE_COMMAND} -E touch ${CUSTOM_BOOST_LOCATION}/boost.built
    WORKING_DIRECTORY ${BOOST_BUILD_DIR}
    DEPENDS ${CUSTOM_BOOST_LOCATION}/boost.configured
    COMMENT "Building Boost")
