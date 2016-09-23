# This ensures that the operator<< overload for Taylor types is available
add_definitions(-DTAYLOR_CXXIO)
#install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/libtaylor DESTINATION ${PROJECT_BINARY_DIR}/external/include)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/libtaylor)
