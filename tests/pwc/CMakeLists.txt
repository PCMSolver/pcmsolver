set(pwc_SOURCES
	pwc.cpp
)

set(test_name pwc)
add_executable(${test_name}.x ${test_name}.cpp )      	
add_dependencies(${test_name}.x googletest)

target_link_libraries(${test_name}.x gtest gtest_main)
target_link_libraries(${test_name}.x solver cavity wavcav green wem pwl utils)
target_link_libraries(${test_name}.x ${LIBS})

add_test(NAME ${test_name} COMMAND ${test_name}.x)

# Export the sources up to src directory level
#set(local_sources "")
#foreach(source ${solver_SOURCES})
#	get_filename_component(sourcepath ${source} ABSOLUTE)
#  	set(local_sources ${sourcepath} ${local_sources})
#endforeach()
#set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)
