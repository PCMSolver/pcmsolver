# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/Config.hpp.in ${PROJECT_BINARY_DIR}/include/Config.hpp @ONLY)

get_property(PCMSOLVER_EXECUTABLE GLOBAL PROPERTY PCMSolver_EXECUTABLE)
# Configure the input parsing script
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in ${PROJECT_BINARY_DIR}/bin/pcmsolver.py @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/bin/pcmsolver.py DESTINATION bin)
