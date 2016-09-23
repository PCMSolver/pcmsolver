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
# autocmake.cfg configuration::
#
# ?

if(NOT DEFINED ${CMAKE_INSTALL_LIBDIR} OR "${${CMAKE_INSTALL_LIBDIR}}" STREQUAL "")
    set(${CMAKE_INSTALL_LIBDIR} lib CACHE STRING "Directory to which libraries installed" FORCE)
endif()

include(GNUInstallDirs)

set(PN ${PROJECT_NAME})

