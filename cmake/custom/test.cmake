option(ENABLE_TESTS "Enable PCMSolver unit tests" ON)

# This macro is used to add a unit test using Boost Unit Testing framework
macro(add_boosttest test_name)
    get_filename_component(the_name ${test_name} NAME_WE)
    add_executable(${the_name}.x ${the_name}.cpp)

    set(_my_libraries "${ARGV1}")
    set(_external_libraries "${ARGV2}")
    # Find threading library aka CMAKE_THREAD_LIBS_INIT
    find_package(Threads)

    # Building on more than one processor can result in race conditions,
    # since custom Boost can be built only on one processor!
    # We thus add this dependency to not get stuck.
    if(BUILD_CUSTOM_BOOST)
        add_dependencies(${the_name}.x custom_boost)
    endif()
    target_link_libraries(${the_name}.x
        ${_my_libraries}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        ${CMAKE_THREAD_LIBS_INIT}
        ${_external_libraries}
        )
    add_test(NAME ${the_name} COMMAND ${the_name}.x)
endmacro()

# This macro is used to copy a file containing reference values into the directory
# where the unit tests will be executed
macro(add_reference reference_file where)
    FILE(COPY ${reference_file} DESTINATION ${where}
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ
        WORLD_READ)
endmacro()

if(ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests) # This must come last!!
endif()
