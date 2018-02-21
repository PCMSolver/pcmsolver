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

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option ${variable}: ${default}")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/tools)
