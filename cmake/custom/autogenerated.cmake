# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/Config.hpp.in
    ${PROJECT_BINARY_DIR}/include/Config.hpp)

get_property(PCMSOLVER_EXECUTABLE GLOBAL PROPERTY PCMSolver_EXECUTABLE)
# Configure the input parsing script
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.py.in pcmsolver.py)
file(COPY ${PROJECT_BINARY_DIR}/pcmsolver.py
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Configure the touch_cmakelists utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/touch_cmakelists.py.in touch_cmakelists.py)
file(COPY ${PROJECT_BINARY_DIR}/touch_cmakelists.py
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Configure the update_gh-pages utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/update_gh-pages.py.in update_gh-pages.py)
file(COPY ${PROJECT_BINARY_DIR}/update_gh-pages.py
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Configure the counter utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/counter.py.in counter.py)
file(COPY ${PROJECT_BINARY_DIR}/counter.py
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Configure the extract_notice utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/extract_notice.py.in extract_notice.py)
file(COPY ${PROJECT_BINARY_DIR}/extract_notice.py
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Configure the update_copyright utility script
configure_file(${PROJECT_SOURCE_DIR}/tools/update_copyright.py.in update_copyright.py)
file(COPY ${PROJECT_BINARY_DIR}/update_copyright.py
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
