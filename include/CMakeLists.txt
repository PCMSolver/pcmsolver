# Configure the header with library-wide preprocessor definitions
configure_file(Config.hpp.in ${PROJECT_SOURCE_DIR}/include/Config.hpp @ONLY)
set_source_files_properties(${PROJECT_SOURCE_DIR}/include/Config.hpp PROPERTIES GENERATED TRUE)
add_custom_target(
  generate-config-hpp
  ALL DEPENDS
    ${PROJECT_SOURCE_DIR}/include/Config.hpp
    ${PROJECT_SOURCE_DIR}/include/Config.hpp.in
    $<$<BOOL:${BUILD_CUSTOM_BOOST}>:custom_boost>
  )

generate_git_info_header(${PROJECT_SOURCE_DIR}/include GitInfo.hpp)

list(APPEND headers_list
  Citation.hpp
  Config.hpp
  Cxx11Workarounds.hpp
  ErrorHandling.hpp
  GitInfo.hpp
  LoggerInterface.hpp
  PhysicalConstants.hpp
  STLUtils.hpp
  TimerInterface.hpp
  )

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
endforeach()
