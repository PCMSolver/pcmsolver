# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/include/Config.hpp.in ${PROJECT_BINARY_DIR}/include/Config.hpp @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/include/Config.hpp DESTINATION include)

# Configure the input parsing script
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in ${PROJECT_BINARY_DIR}/bin/pcmsolver.py @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/bin/pcmsolver.py DESTINATION bin)
# Install GetKw Python bindings
install(FILES ${PROJECT_SOURCE_DIR}/tools/getkw.py
              ${PROJECT_SOURCE_DIR}/tools/pyparsing.py
        DESTINATION bin)

# Install docopt.py in the bin subdirectory
file(COPY ${PROJECT_SOURCE_DIR}/cmake/lib/docopt/docopt.py DESTINATION ${PROJECT_BINARY_DIR}/bin)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/lib/docopt/docopt.py DESTINATION bin)
