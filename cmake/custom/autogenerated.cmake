if(PYMOD_INSTALL_LIBDIR)
    set(PYMOD_INSTALL_FULLDIR "${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pcmsolver")
else()
    set(PYMOD_INSTALL_FULLDIR "${CMAKE_INSTALL_BINDIR}")
endif()

# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/include/Config.hpp.in
               ${PROJECT_BINARY_DIR}/include/tmp-config-hpp @ONLY)
add_custom_command(
  DEPENDS ${PROJECT_SOURCE_DIR}/include/Config.hpp.in
  OUTPUT
    ${PROJECT_BINARY_DIR}/include/Config.hpp
  COMMAND
    cmake -E copy ${PROJECT_BINARY_DIR}/include/tmp-config-hpp ${PROJECT_BINARY_DIR}/include/Config.hpp
  VERBATIM
  )
set_source_files_properties(${PROJECT_BINARY_DIR}/include/Config.hpp PROPERTIES GENERATED TRUE)
add_custom_target(generate-config-hpp ALL DEPENDS ${PROJECT_BINARY_DIR}/include/Config.hpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(generate-config-hpp custom_boost)
endif()
install(FILES ${PROJECT_BINARY_DIR}/include/Config.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

generate_git_info_header(${PROJECT_BINARY_DIR}/include GitInfo.hpp)
install(FILES ${PROJECT_BINARY_DIR}/include/GitInfo.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

# Configure the input parsing script
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/tmp-pcmsolver-py @ONLY)
add_custom_command(
  DEPENDS ${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in
  OUTPUT
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py
  COMMAND
    cmake -E copy ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/tmp-pcmsolver-py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py
  VERBATIM
  )
add_custom_target(generate-pcmsolver-py ALL DEPENDS ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py)
install(FILES ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmsolver.py DESTINATION ${PYMOD_INSTALL_FULLDIR})
# Configure the codata Python module
configure_file(${PROJECT_SOURCE_DIR}/tools/codata.py.in ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/codata.py @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/codata.py DESTINATION ${PYMOD_INSTALL_FULLDIR})

# Install GetKw Python bindings
# If using Python 3 use py3.x-getkw.py
if(PYTHON_VERSION_MAJOR VERSION_EQUAL 3)
  configure_file(${PROJECT_SOURCE_DIR}/tools/py3.x-getkw.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py COPYONLY)
else()
  configure_file(${PROJECT_SOURCE_DIR}/tools/py2.x-getkw.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py COPYONLY)
endif()
install(FILES ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py
  DESTINATION ${PYMOD_INSTALL_FULLDIR})
file(COPY ${PROJECT_SOURCE_DIR}/tools/pyparsing.py
     DESTINATION ${PYMOD_INSTALL_FULLDIR})
install(FILES ${PROJECT_SOURCE_DIR}/tools/pyparsing.py
        DESTINATION ${PYMOD_INSTALL_FULLDIR})

# Install docopt.py in the bin subdirectory
file(COPY ${PROJECT_SOURCE_DIR}/cmake/lib/docopt/docopt.py
     DESTINATION ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})
install(FILES ${PROJECT_SOURCE_DIR}/cmake/lib/docopt/docopt.py
        DESTINATION ${PYMOD_INSTALL_FULLDIR})
