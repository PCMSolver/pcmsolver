# Configure the header with library-wide preprocessor definitions
configure_file(Config.hpp.in ${PROJECT_BINARY_DIR}/include/Config.hpp @ONLY)
set_source_files_properties(${PROJECT_BINARY_DIR}/include/Config.hpp PROPERTIES GENERATED TRUE)
add_custom_target(
  generate-config-hpp
  ALL DEPENDS
    ${PROJECT_BINARY_DIR}/include/Config.hpp
    ${PROJECT_SOURCE_DIR}/include/Config.hpp.in
    $<$<BOOL:${BUILD_CUSTOM_BOOST}>:custom_boost>
  )

install(
  FILES
    Citation.hpp
    ${PROJECT_BINARY_DIR}/include/Config.hpp
    ErrorHandling.hpp
    LoggerInterface.hpp
    PhysicalConstants.hpp
    STLUtils.hpp
    TimerInterface.hpp
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  )
