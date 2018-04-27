configure_file(go_pcm.py.in ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py @ONLY)
file(
  COPY
    plot_cavity.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/plot_cavity.py
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

configure_file(init.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py COPYONLY)
file(
  COPY
    getkw.py
    pcmdata.py
    pcmparser.py
    pyparsing.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmdata.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmparser.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pyparsing.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )

# Install docopt.py provided by Autocmake
file(TO_NATIVE_PATH "${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/docopt" PYMOD_INSTALL_EXTDIR)
file(
  COPY
    ${PROJECT_SOURCE_DIR}/cmake/autocmake/external/__init__.py
    ${PROJECT_SOURCE_DIR}/cmake/autocmake/external/docopt.py
  DESTINATION
  ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/external
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/external/__init__.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/external/docopt.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}/external
  )

