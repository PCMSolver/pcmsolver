cmake_policy(SET CMP0022 NEW)
cmake_policy(SET CMP0063 NEW)

include_directories(
  ${CMAKE_CURRENT_LIST_DIR}
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include
  SYSTEM ${CMAKE_CURRENT_LIST_DIR}/utils/getkw
  )

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(bi_operators)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)

list(APPEND _objects
    $<TARGET_OBJECTS:cavity>
    $<$<BOOL:${ENABLE_Fortran_API}>:$<TARGET_OBJECTS:fortran_bindings>>
    $<TARGET_OBJECTS:green>
    $<TARGET_OBJECTS:interface>
    $<TARGET_OBJECTS:metal>
    $<TARGET_OBJECTS:bi_operators>
    $<TARGET_OBJECTS:pedra>
    $<TARGET_OBJECTS:solver>
    $<TARGET_OBJECTS:utils>
    $<TARGET_OBJECTS:getkw>
    )

if(NOT STATIC_LIBRARY_ONLY)
  add_library(pcm-shared SHARED ${_objects})
  target_compile_definitions(pcm-shared
    PUBLIC
      PCMSolver_EXPORTS
    INTERFACE
      $<INSTALL_INTERFACE:USING_PCMSolver>
    )
  target_include_directories(pcm-shared
    INTERFACE
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  target_link_libraries(pcm-shared
    PRIVATE
      ${LIBC_INTERJECT}
      ZLIB::ZLIB
    )
  set_target_properties(pcm-shared
    PROPERTIES
      SOVERSION 1
      CXX_VISIBILITY_PRESET hidden
      VISIBILITY_INLINES_HIDDEN 1
      MACOSX_RPATH ON
      OUTPUT_NAME "pcm"
      EXPORT_NAME "pcm"
    )
  install(
    TARGETS
      pcm-shared
    EXPORT
      "PCMSolverTargets-shared"
    RUNTIME DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    )
  add_dependencies(pcm-shared generate-config-hpp)
endif()

if(NOT SHARED_LIBRARY_ONLY)
  add_library(pcm-static STATIC ${_objects})
  target_compile_definitions(pcm-static
    PUBLIC
      PCMSolver_STATIC_DEFINE
    INTERFACE
      $<INSTALL_INTERFACE:USING_PCMSolver>
    )
  target_include_directories(pcm-static
    INTERFACE
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  target_link_libraries(pcm-static
    PRIVATE
      ZLIB::ZLIB
    )
  set_target_properties(pcm-static
    PROPERTIES
       OUTPUT_NAME "pcm"
       EXPORT_NAME "pcm"
    )
  install(
    TARGETS
      pcm-static
    EXPORT
      "PCMSolverTargets-static"
    ARCHIVE DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    )
  add_dependencies(pcm-static generate-config-hpp)
endif()

if(STATIC_LIBRARY_ONLY)
  add_library(PCMSolver ALIAS pcm-static)
else()
  add_library(PCMSolver ALIAS pcm-shared)
endif()

include(GenerateExportHeader)
generate_export_header(PCMSolver
  BASE_NAME "PCMSolver"
  EXPORT_MACRO_NAME "PCMSolver_EXPORT"
  EXPORT_FILE_NAME "PCMSolverExport.h"
  DEPRECATED_MACRO_NAME "PCMSolver_DEPRECATED"
  NO_EXPORT_MACRO_NAME "PCMSolver_NO_EXPORT"
  STATIC_DEFINE "PCMSolver_STATIC_DEFINE"
  NO_DEPRECATED_MACRO_NAME "PCMSolver_NO_DEPRECATED"
  DEFINE_NO_DEPRECATED)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/PCMSolverExport.h DESTINATION ${PROJECT_BINARY_DIR}/include)
install(FILES ${PROJECT_BINARY_DIR}/include/PCMSolverExport.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PCMSolver)

# This directory contains source file for executables
if(BUILD_STANDALONE)
  add_subdirectory(bin)
endif()

# <<<  Export Config  >>>

include(CMakePackageConfigHelpers)
# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")
configure_package_config_file(
  ${PCMSolver_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

add_custom_target(update_version
  ALL
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/versioner.py
        --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        --cmakeout ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
  COMMAND
    ${CMAKE_COMMAND}
         -DWTO="${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}"
         -DPN="PCMSolver"
         -P ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT
    "Generating version info"
  )

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )
if(NOT STATIC_LIBRARY_ONLY)
  install(
    EXPORT
      "PCMSolverTargets-shared"
    NAMESPACE
      "PCMSolver::"
    DESTINATION
      ${CMAKECONFIG_INSTALL_DIR}
    )
endif()

if(NOT SHARED_LIBRARY_ONLY)
  install(
    EXPORT
      "PCMSolverTargets-static"
    NAMESPACE
      "PCMSolver::"
    DESTINATION
      ${CMAKECONFIG_INSTALL_DIR}
    )
endif()
