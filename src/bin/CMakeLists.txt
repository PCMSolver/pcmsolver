# Executables for testing separately from the unit test framework or for profiling purposes

#add_executable(debug_wavcav debug_wavcav.cpp)
#target_link_libraries(debug_wavcav solver cavity wavcav pwl wem utils)

add_executable(run_pcm run_pcm.cpp)
add_executable(PCMSolver::run_pcm ALIAS run_pcm)
target_link_libraries(run_pcm PCMSolver::pcm)
if(ENABLE_GPERFTOOLS)
  target_link_libraries(run_pcm ${GPERFTOOLS_LIBRARIES})
endif()
target_compile_definitions(run_pcm
  PRIVATE
    "$<$<BOOL:${HAS_CXX11}>:HAS_CXX11>"
  )
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
    EXPORT_NAME "exe"
  )
install(
  TARGETS
    run_pcm
  EXPORT
    executable_interface
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )
install(
  EXPORT
    executable_interface
  FILE
    "PCMSolverTargets-exe"
  NAMESPACE
    "PCMSolver::"
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )
