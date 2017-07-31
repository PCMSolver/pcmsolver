# Executables for testing separate from the unit test framework or for profiling purposes
# The executables will be created in ${PROJECT_BINARY_DIR}/src/bin

#add_executable(debug_wavcav debug_wavcav.cpp)
#target_link_libraries(debug_wavcav solver cavity wavcav pwl wem utils)

link_directories(${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
add_executable(run_pcm run_pcm.cpp)
add_dependencies(run_pcm generate-config-hpp)
if(STATIC_LIBRARY_ONLY)
  target_link_libraries(run_pcm pcm-static)
else()
  target_link_libraries(run_pcm pcm-shared)
endif()
target_compile_options(run_pcm PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")
set_target_properties(run_pcm PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set_target_properties(run_pcm PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
set_target_properties(run_pcm PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR};${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
install(TARGETS run_pcm
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
