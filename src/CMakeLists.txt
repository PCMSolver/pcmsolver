cmake_policy (SET CMP0022 NEW)

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
    target_link_libraries(pcm-shared PRIVATE ${LIBC_INTERJECT}
                                             ZLIB::ZLIB)
    set_target_properties(pcm-shared PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR}
                                                CXX_VISIBILITY_PRESET hidden
                                                VISIBILITY_INLINES_HIDDEN 1
                                                OUTPUT_NAME "pcm"
                                                EXPORT_NAME "pcm")

    install(TARGETS pcm-shared
            EXPORT "${PROJECT_NAME}Targets-shared"
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    add_dependencies(pcm-shared generate-config-hpp)
endif()

if(NOT SHARED_LIBRARY_ONLY)
    add_library(pcm-static STATIC ${_objects})
    target_link_libraries(pcm-static PRIVATE ZLIB::ZLIB)
    set_target_properties(pcm-static PROPERTIES COMPILE_FLAGS -D${PROJECT_NAME}_STATIC_DEFINE
                                                OUTPUT_NAME "pcm"
                                                EXPORT_NAME "pcm")
    install(TARGETS pcm-static
            EXPORT "${PROJECT_NAME}Targets-static"
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    add_dependencies(pcm-static generate-config-hpp)
endif()

if(STATIC_LIBRARY_ONLY)
    add_library(${PROJECT_NAME} ALIAS pcm-static)
else ()
    add_library(${PROJECT_NAME} ALIAS pcm-shared)
endif()

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

# <<<  Export interface  >>>

if(NOT STATIC_LIBRARY_ONLY)
    if (APPLE)
        set_target_properties(pcm-shared PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()
    target_compile_definitions(pcm-shared INTERFACE USING_${PROJECT_NAME})
    target_include_directories(pcm-shared INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
endif()

if(NOT SHARED_LIBRARY_ONLY)
    target_compile_definitions(pcm-static INTERFACE USING_${PROJECT_NAME})
    target_include_directories(pcm-static INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
endif()

# <<<  Export Config  >>>

    # explicit "share" not "DATADIR" for CMake search path
set (CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")
if(NOT STATIC_LIBRARY_ONLY)
    install (EXPORT "${PROJECT_NAME}Targets-shared"
             NAMESPACE "${PROJECT_NAME}::"
             DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()

if(NOT SHARED_LIBRARY_ONLY)
    install (EXPORT "${PROJECT_NAME}Targets-static"
             NAMESPACE "${PROJECT_NAME}::"
             DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()
