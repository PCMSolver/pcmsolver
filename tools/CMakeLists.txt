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

# Install the codata Python module
file(COPY ${PROJECT_SOURCE_DIR}/tools/codata.py DESTINATION ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})
install(FILES ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/codata.py DESTINATION ${PYMOD_INSTALL_FULLDIR})

# Install the plot_cavity Python script
file(COPY ${PROJECT_SOURCE_DIR}/tools/plot_cavity.py DESTINATION ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})
install(FILES ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/plot_cavity.py DESTINATION ${PYMOD_INSTALL_FULLDIR})

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

# Install __init__.py
configure_file(${PROJECT_SOURCE_DIR}/tools/init.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py COPYONLY)
install(FILES ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
        DESTINATION ${PYMOD_INSTALL_FULLDIR})
