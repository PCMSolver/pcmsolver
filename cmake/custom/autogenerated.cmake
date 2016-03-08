# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/include/Config.hpp.in
               ${PROJECT_BINARY_DIR}/include/tmp-config-hpp @ONLY)
add_custom_command(
  OUTPUT
    ${PROJECT_BINARY_DIR}/include/Config.hpp
  COMMAND
    cmake -E copy ${PROJECT_BINARY_DIR}/include/tmp-config-hpp ${PROJECT_BINARY_DIR}/include/Config.hpp
  VERBATIM
  )
set_source_files_properties(${PROJECT_BINARY_DIR}/include/Config.hpp PROPERTIES GENERATED TRUE)
add_custom_target(generate-config-hpp ALL DEPENDS ${PROJECT_BINARY_DIR}/include/Config.hpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(generate-config-hpp custom_boost)
endif()
install(FILES ${PROJECT_BINARY_DIR}/include/Config.hpp DESTINATION include)


# Configure the input parsing script
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in ${PROJECT_BINARY_DIR}/bin/pcmsolver.py @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/bin/pcmsolver.py DESTINATION bin)
# Install GetKw Python bindings
file(COPY ${PROJECT_SOURCE_DIR}/tools/getkw.py
          ${PROJECT_SOURCE_DIR}/tools/pyparsing.py
     DESTINATION bin)
install(FILES ${PROJECT_SOURCE_DIR}/tools/getkw.py
              ${PROJECT_SOURCE_DIR}/tools/pyparsing.py
        DESTINATION bin)

# Install docopt.py in the bin subdirectory
file(COPY ${PROJECT_SOURCE_DIR}/cmake/lib/docopt/docopt.py DESTINATION ${PROJECT_BINARY_DIR}/bin)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/lib/docopt/docopt.py DESTINATION bin)
