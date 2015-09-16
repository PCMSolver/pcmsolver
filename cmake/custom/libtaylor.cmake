set(TaylorCMakeArgs
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
   )
# This ensures that the operator<< overload for Taylor types is available
add_definitions(-DTAYLOR_CXXIO)

ExternalProject_Add(libtaylor
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${TaylorCMakeArgs}
    )
