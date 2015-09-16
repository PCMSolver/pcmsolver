include(ExternalProject)

set(GetkwCMakeArgs
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
   -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
   -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}
   )

ExternalProject_Add(libgetkw
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libgetkw
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libgetkw-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/libgetkw-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/libgetkw-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${GetkwCMakeArgs}
    )

if(BUILD_CUSTOM_BOOST)
	add_dependencies(libgetkw custom_boost)
endif()

link_directories(${PROJECT_BINARY_DIR}/external/lib)

set(GETKW_PYTHON_DIR ${PROJECT_BINARY_DIR}/external/share/libgetkw)
