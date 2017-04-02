# This ensures that the operator<< overload for Taylor types is available
add_definitions(-DTAYLOR_CXXIO)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/libtaylor
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/external)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/libtaylor)
