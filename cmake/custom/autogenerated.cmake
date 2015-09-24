# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/Config.hpp.in ${PROJECT_BINARY_DIR}/include/Config.hpp @ONLY)

get_property(PCMSOLVER_EXECUTABLE GLOBAL PROPERTY PCMSolver_EXECUTABLE)
# Configure the input parsing script
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in ${PROJECT_BINARY_DIR}/bin/pcmsolver.py @ONLY)

# Configure the extract_notice utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/extract_notice.py.in ${PROJECT_BINARY_DIR}/bin/extract_notice.py @ONLY)

# Configure the update_copyright utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/update_copyright.py.in ${PROJECT_BINARY_DIR}/bin/update_copyright.py @ONLY)
