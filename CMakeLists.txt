cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(PCMSolver CXX C Fortran)

option(BUILD_AS_LIBRARY   "Build as library"                                    OFF)
option(ENABLE_EIGEN_MKL   "Enable Eigen3 automatic fallback to some of Intel MKL algorithms"                 OFF)

set(PCMSolver_VERSION_MAJOR 0)
set(PCMSolver_VERSION_MINOR 1)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    )

set(LIBS)

include(ConfigVersion)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
#include(ConfigDocumentation)
include(ConfigExternal)
include(ConfigMath)
include(ConfigSafeGuards)
include(GenericMacros)

find_package(Eigen3 3.0.0)

add_custom_target(external)

set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DPARENT_INCLUDE_DIR=${PROJECT_BINARY_DIR}/external/include
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
        -DPARENT_BINARY_DIR=${PROJECT_BINARY_DIR}
   )

add_external(libgetkw)

ExternalProject_Add(libtaylor
	SVN_REPOSITORY http://libtaylor.googlecode.com/svn/trunk/
       	DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
	)

add_dependencies(external libgetkw libtaylor)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR}/external/libtaylor-build)
include_directories(${PROJECT_BINARY_DIR}/external/libtaylor-build/modules)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external/include)

link_directories(${PROJECT_BINARY_DIR}/external/lib)
link_directories(${PROJECT_BINARY_DIR}/external/libtaylor-build/external/lib)

add_subdirectory(src)

include_directories (src/solverlib)
include_directories (src/wavcavlib)
include_directories (src/wemlib)
include_directories (src/pwllib)
include_directories (src/pedralib)
include_directories (src/cavitylib)
include_directories (src/metal)
include_directories (src/green)

add_executable(wavcav.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wavcav.cpp)
target_link_libraries(wavcav.x cavity pedra wavcav pwl wem getkw ${LIBS})

add_executable(pcmsolver.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/pcmsolver.cpp)
target_link_libraries(pcmsolver.x solver green cavity pedra metal pwl wem wavcav getkw ${LIBS})

add_executable(pwl.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/pwllib/main_pwl.c)
target_link_libraries(pwl.x pwl wem ${LIBS})

add_executable(pwc.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wemlib/main.c)
target_link_libraries(pwc.x wem ${LIBS})

add_executable(cav.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wavcavlib/main.c)
target_link_libraries(cav.x wavcav ${LIBS})

