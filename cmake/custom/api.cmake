#.rst:
#
# Manage compilation of API.
# Optionally, enable compilation of Fortran 90 API bindings.
#
# Variables defined::
#
#   ENABLE_FORTRAN_API
#
# autocmake.yml configuration::
#
#   docopt: "--fbindings=<ENABLE_FORTRAN_API> Enable compilation of Fortran 90 API bindings <ON/OFF> [default: ON]."
#   define: "'-DENABLE_FORTRAN_API={0}'.format(arguments['--fbindings'])"

option(ENABLE_FORTRAN_API "Enable compilation of Fortran 90 API bindings" ON)

add_subdirectory(api)
include_directories(${PROJECT_SOURCE_DIR}/api)
