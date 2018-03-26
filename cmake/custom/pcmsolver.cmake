option(ENABLE_LOGGER "Enable logger" OFF)
option(ENABLE_TIMER "Enable timer" ON)
option(BUILD_STANDALONE "Enable build of standalone executables" ON)

# Add definitions
if(ENABLE_TIMER)
  add_definitions(-DENABLE_TIMER)
endif()
if(ENABLE_LOGGER)
  add_definitions(-DENABLE_LOGGER)
endif()

set(BOOST_MINIMUM_REQUIRED 1.54.0)
set(BOOST_COMPONENTS_REQUIRED "")

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/pcmsolver" PYMOD_INSTALL_FULLDIR)

add_custom_target(update_version
  ALL
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/versioner.py
        --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        --cmakeout ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
        --headerout ${CMAKE_CURRENT_BINARY_DIR}/include/VersionInfo.hpp
  COMMAND
    ${CMAKE_COMMAND}
         -DWTO="${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}"
         -DPN="PCMSolver"
         -P ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT
    "Generating version info"
  )
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/include/VersionInfo.hpp
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  )

add_subdirectory(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/tools)
