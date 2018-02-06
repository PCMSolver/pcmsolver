if(PYMOD_INSTALL_LIBDIR)
  set(PYMOD_INSTALL_FULLDIR "${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pcmsolver")
else()
  set(PYMOD_INSTALL_FULLDIR "${CMAKE_INSTALL_LIBDIR}/python/pcmsolver")
endif()

add_subdirectory(tools)
