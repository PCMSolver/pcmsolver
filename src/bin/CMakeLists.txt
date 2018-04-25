# Executables for testing separately from the unit test framework or for profiling purposes

#add_executable(debug_wavcav debug_wavcav.cpp)
#target_link_libraries(debug_wavcav solver cavity wavcav pwl wem utils)

add_executable(run_pcm run_pcm.cpp)
target_link_libraries(run_pcm PCMSolver)
if(ENABLE_GPERFTOOLS)
  target_link_libraries(run_pcm ${GPERFTOOLS_LIBRARIES})
endif()
target_compile_options(run_pcm
  PRIVATE
    "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>"
  )
set_target_properties(run_pcm
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH ON
  )
install(TARGETS run_pcm RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
