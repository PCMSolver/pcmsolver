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
#   docopt:
#     - "--bindir=<CMAKE_INSTALL_BINDIR> User executables [default: bin]."
#     - "--libdir=<CMAKE_INSTALL_LIBDIR> Object code libraries [default: lib]."
#     - "--includedir=<CMAKE_INSTALL_INCLUDEDIR> C header files [default: include]."
#     - "--datadir=<CMAKE_INSTALL_DATADIR> Read-only architecture-independent data root [default: share]."
#   define:
#     - "'-DCMAKE_INSTALL_BINDIR={0}'.format(arguments['--bindir'])"
#     - "'-DCMAKE_INSTALL_LIBDIR={0}'.format(arguments['--libdir'])"
#     - "'-DCMAKE_INSTALL_INCLUDEDIR={0}'.format(arguments['--includedir'])"
#     - "'-DCMAKE_INSTALL_DATADIR={0}'.format(arguments['--datadir'])"

include(GNUInstallDirs)
