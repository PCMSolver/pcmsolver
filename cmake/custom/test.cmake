option(ENABLE_TESTS "Enable PCMSolver unit tests" ON)

if(ENABLE_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(tests) # This must come last!!
endif()
