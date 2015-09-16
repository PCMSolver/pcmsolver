option(ENABLE_LOGGER "Enable logger" ON)
option(ENABLE_TIMER "Enable timer" ON)
option(BUILD_STANDALONE "Enable build of standalone executables" ON)
option(ENABLE_FORTRAN_API "Builds optional Fortran90 API" OFF)

# PCMSolver sub-libraries
set_property(GLOBAL PROPERTY PCMSolver_LIBRARIES)
# PCMSolver headers
set_property(GLOBAL PROPERTY PCMSolver_HEADER_DIRS)

include_directories(SYSTEM ${PROJECT_BINARY_DIR}/external/include)
include_directories(${PROJECT_BINARY_DIR}/include)
