include_directories(${CMAKE_CURRENT_LIST_DIR}/cavity
                    ${CMAKE_CURRENT_LIST_DIR}/green
                    ${CMAKE_CURRENT_LIST_DIR}/green/dielectric_profile
                    ${CMAKE_CURRENT_LIST_DIR}/interface
                    ${CMAKE_CURRENT_LIST_DIR}/bi_operators
                    ${CMAKE_CURRENT_LIST_DIR}/solver
                    ${CMAKE_CURRENT_LIST_DIR}/utils
                   )

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(bi_operators)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)

get_property(_pcmsolver_cxx_sources GLOBAL PROPERTY PCMSolver_CXX_SOURCES)
get_property(_pcmsolver_c_sources GLOBAL PROPERTY PCMSolver_C_SOURCES)
get_property(_pcmsolver_fortran_sources GLOBAL PROPERTY PCMSolver_Fortran_SOURCES)
add_library(pcm-base OBJECT ${_pcmsolver_cxx_sources} ${_pcmsolver_c_sources} ${_pcmsolver_fortran_sources})
add_dependencies(pcm-base libgetkw)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(pcm-base custom_boost)
endif()
if(NOT STATIC_LIBRARY_ONLY)
  set_target_properties(pcm-base PROPERTIES POSITION_INDEPENDENT_CODE 1
                                            INTERPROCEDURAL_OPTIMIZATION 1)
  add_library(pcm-shared SHARED $<TARGET_OBJECTS:pcm-base>)
  target_link_libraries(pcm-shared ${SUBMODULES_INSTALL_PREFIX}/lib/libgetkw.a ${ZLIB_LIBRARIES})
  set_target_properties(pcm-shared PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME "pcm")
  install(TARGETS pcm-shared DESTINATION lib)
endif()

add_library(pcm-static STATIC $<TARGET_OBJECTS:pcm-base>)
set_target_properties(pcm-static PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME "pcm")
install(TARGETS pcm-static ARCHIVE DESTINATION lib)

# This directory contains source file for executables
if(BUILD_STANDALONE)
   add_subdirectory(bin)
endif()
