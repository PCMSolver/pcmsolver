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

include_directories(${PROJECT_SOURCE_DIR}/api)

install(FILES ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(FILES ${PROJECT_SOURCE_DIR}/api/PCMInput.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

if(ENABLE_FORTRAN_API)
    add_definitions(-DENABLE_FORTRAN_API)
    add_library(fortran_bindings OBJECT ${PROJECT_SOURCE_DIR}/api/pcmsolver.f90)
    list(APPEND _objects $<TARGET_OBJECTS:fortran_bindings>)
endif()
