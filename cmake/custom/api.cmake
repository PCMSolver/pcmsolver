#.rst:
#
# Manage compilation of API.
# Optionally, enable testing of Fortran 90 API bindings.
#
# Variables defined::
#
#   TEST_Fortran_API
#
# autocmake.yml configuration::
#
#   docopt: "--fbindings=<TEST_Fortran_API> Enable testing of Fortran 90 API bindings <ON/OFF> [default: ON]."
#   define: "'-DTEST_Fortran_API={0}'.format(arguments['--fbindings'])"

option(TEST_Fortran_API "Enable testing of Fortran 90 API bindings" ON)

add_subdirectory(api)
include_directories(${PROJECT_SOURCE_DIR}/api)
