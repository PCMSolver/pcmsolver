option(ENABLE_LOGGER "Enable logger" ON)
option(ENABLE_TIMER "Enable timer" ON)
option(BUILD_STANDALONE "Enable build of standalone executables" ON)
option(ENABLE_FORTRAN_API "Builds optional Fortran90 API" OFF)

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

include_directories(SYSTEM ${PROJECT_BINARY_DIR}/external/include)
include_directories(${PROJECT_BINARY_DIR}/include)
