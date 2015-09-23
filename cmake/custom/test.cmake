option(ENABLE_TESTS "Enable PCMSolver unit tests" ON)

# This macro is used to copy a file containing reference values into the directory
# where the unit tests will be executed
macro(add_reference reference_file where)
    FILE(COPY ${reference_file} DESTINATION ${where}
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ
        WORLD_READ)
endmacro()

if(ENABLE_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(tests) # This must come last!!
endif()
