set(wavelet_cavity_SOURCES
	wavelet_cavity.cpp
)

set(test_name wavelet_cavity)
add_executable(${test_name}.x ${test_name}.cpp)      	
add_dependencies(${test_name}.x googletest)

target_link_libraries(${test_name}.x gtest gtest_main pthread)
target_link_libraries(${test_name}.x solver cavity wavcav green wem pwl utils)
target_link_libraries(${test_name}.x ${LIBS})

#target_link_libraries(${test_name}.x
#			solver cavity green pedra tsless utils
#                      	${LIBS}
#                      	${GTEST_LIBS_DIR}/libgtest.a
#                      	${GTEST_LIBS_DIR}/libgtest_main.a
#                     )
add_test(NAME ${test_name} COMMAND ${test_name}.x)
#add_library(solver ${solver_SOURCES})
#add_dependencies(solver cavity green wem pwl)

# Export the sources up to src directory level
#set(local_sources "")
#foreach(source ${solver_SOURCES})
#	get_filename_component(sourcepath ${source} ABSOLUTE)
#  	set(local_sources ${sourcepath} ${local_sources})
#endforeach()
#set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)
