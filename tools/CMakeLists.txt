file(
  COPY
    go_pcm.py
    plot_cavity.py
    setup.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/plot_cavity.py
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/setup.py
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

if(PYTHON_VERSION_MAJOR VERSION_EQUAL 3)
  configure_file(py3.x-getkw.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py COPYONLY)
else()
  configure_file(py2.x-getkw.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py COPYONLY)
endif()
configure_file(init.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py COPYONLY)
file(
  COPY
    codata.py
    pcmparser.py
    pyparsing.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/codata.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmparser.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pyparsing.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
