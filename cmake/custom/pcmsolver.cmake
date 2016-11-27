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

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/src/utils/getkw)
