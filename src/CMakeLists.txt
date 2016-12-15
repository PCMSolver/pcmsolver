include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories(SYSTEM ${CMAKE_CURRENT_LIST_DIR}/utils/getkw)

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(bi_operators)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)

list(APPEND _objects $<TARGET_OBJECTS:cavity>
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
  add_library(${PROJECT_NAME} ALIAS pcm-shared)
  target_link_libraries(pcm-shared ${ZLIB_LIBRARIES} ${LIBC_INTERJECT})
  set_target_properties(pcm-shared PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR}
                                              CXX_VISIBILITY_PRESET hidden
                                              VISIBILITY_INLINES_HIDDEN 1
                                              OUTPUT_NAME "pcm")
  install(TARGETS pcm-shared
          RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  add_dependencies(pcm-shared generate-config-hpp)
endif()

add_library(pcm-static STATIC ${_objects})
if(NOT TARGET ${PROJECT_NAME})
  add_library(${PROJECT_NAME} ALIAS pcm-static)
endif()
set_target_properties(pcm-static PROPERTIES COMPILE_FLAGS -D${PROJECT_NAME}_STATIC_DEFINE
                                            OUTPUT_NAME "pcm")
install(TARGETS pcm-static
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
add_dependencies(pcm-static generate-config-hpp)

include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME}
  BASE_NAME "PCMSolver"
  EXPORT_MACRO_NAME "PCMSolver_EXPORT"
  EXPORT_FILE_NAME "PCMSolverExport.h"
  DEPRECATED_MACRO_NAME "PCMSolver_DEPRECATED"
  NO_EXPORT_MACRO_NAME "PCMSolver_NO_EXPORT"
  STATIC_DEFINE "PCMSolver_STATIC_DEFINE"
  NO_DEPRECATED_MACRO_NAME "PCMSolver_NO_DEPRECATED"
  DEFINE_NO_DEPRECATED)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/PCMSolverExport.h DESTINATION ${PROJECT_BINARY_DIR}/include)
install(FILES ${PROJECT_BINARY_DIR}/include/PCMSolverExport.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

# This directory contains source file for executables
if(BUILD_STANDALONE)
   add_subdirectory(bin)
endif()
