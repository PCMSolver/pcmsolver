option(ENABLE_LOGGER "Enable logger" OFF)
option(ENABLE_TIMER "Enable timer" ON)
option(BUILD_STANDALONE "Enable build of standalone executables" ON)
option(ENABLE_FORTRAN_API "Builds optional Fortran90 API" OFF)

# Add definitions
if(ENABLE_TIMER)
  add_definitions(-DENABLE_TIMER)
endif()
if(ENABLE_LOGGER)
  add_definitions(-DENABLE_LOGGER)
endif()

# This can be set by the host project
# and tweaks the location of the submodules install location
if(NOT DEFINED SUBMODULES_INSTALL_PREFIX)
    set(SUBMODULES_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
endif()

set(BOOST_MINIMUM_REQUIRED 1.54.0)
set(BOOST_COMPONENTS_REQUIRED "")
if(BUILD_STANDALONE OR ENABLE_TIMER)
  set(BOOST_COMPONENTS_REQUIRED "chrono;timer;system")
endif()

# librt : needed for the timers (Linux only)
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  find_package(RT REQUIRED)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# PCMSolver C++ sources
set_property(GLOBAL PROPERTY PCMSolver_CXX_SOURCES)
# PCMSolver C sources
set_property(GLOBAL PROPERTY PCMSolver_C_SOURCES)
# PCMSolver Fortran sources
set_property(GLOBAL PROPERTY PCMSolver_Fortran_SOURCES)
# PCMSolver headers
set_property(GLOBAL PROPERTY PCMSolver_HEADER_DIRS)

include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(SYSTEM ${SUBMODULES_INSTALL_PREFIX}/include)
