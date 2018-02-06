# Configure the header with library-wide preprocessor definitions
configure_file(${PROJECT_SOURCE_DIR}/include/Config.hpp.in
               ${PROJECT_BINARY_DIR}/include/tmp-config-hpp @ONLY)
add_custom_command(
  DEPENDS ${PROJECT_SOURCE_DIR}/include/Config.hpp.in
  OUTPUT
    ${PROJECT_BINARY_DIR}/include/Config.hpp
  COMMAND
    cmake -E copy ${PROJECT_BINARY_DIR}/include/tmp-config-hpp ${PROJECT_BINARY_DIR}/include/Config.hpp
  VERBATIM
  )
set_source_files_properties(${PROJECT_BINARY_DIR}/include/Config.hpp PROPERTIES GENERATED TRUE)
add_custom_target(generate-config-hpp ALL DEPENDS ${PROJECT_BINARY_DIR}/include/Config.hpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(generate-config-hpp custom_boost)
endif()
install(FILES ${PROJECT_BINARY_DIR}/include/Config.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

generate_git_info_header(${PROJECT_BINARY_DIR}/include GitInfo.hpp)
install(FILES ${PROJECT_BINARY_DIR}/include/GitInfo.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

list(APPEND headers_list Citation.hpp
                         Cxx11Workarounds.hpp
                         ErrorHandling.hpp
                         LoggerInterface.hpp
                         PhysicalConstants.hpp
                         STLUtils.hpp
                         TimerInterface.hpp
                         )

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
endforeach()
