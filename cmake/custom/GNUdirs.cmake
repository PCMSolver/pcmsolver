#.rst:
#
# Enables manipulation of install directories.
#
# Variables modified::
#
#   CMAKE_INSTALL_BINDIR
#   CMAKE_INSTALL_LIBDIR
#   CMAKE_INSTALL_INCLUDEDIR
#   CMAKE_INSTALL_DATADIR
#
# autocmake.yml configuration::
#
# ?

if("${${CMAKE_INSTALL_LIBDIR}}" STREQUAL "" OR NOT DEFINED ${CMAKE_INSTALL_LIBDIR})
  set(${CMAKE_INSTALL_LIBDIR} lib CACHE STRING "Directory to which libraries are installed" FORCE)
endif()

include(GNUInstallDirs)
