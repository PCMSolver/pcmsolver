cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(PCMSolver)
enable_language(CXX C Fortran)

option(ENABLE_EIGEN_MKL   "Enable Eigen3 automatic fallback to some of Intel MKL algorithms"                 OFF)

set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 1)
set (PROJECT_VERSION_PATCH 0)

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
include(ConfigDocumentation)
include(ConfigExternal)
include(ConfigMath)
include(ConfigSafeGuards)
include(GenericMacros)
include(ConfigGitRevision)


find_package(Eigen3 3.0.0)
if(ENABLE_EIGEN_MKL)
	message("-- This option requires at least Eigen 3.1.0 and Intel MKL 10.3, be sure you have read http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html")
        set(EIGEN_USE_MKL_ALL ON)
endif()

find_package(PythonInterp REQUIRED)

add_custom_target(external)

set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE="Release"
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

configure_file (
    ${PROJECT_SOURCE_DIR}/Config.h.in
    ${PROJECT_SOURCE_DIR}/Config.h
    )

install(FILES Config.h DESTINATION include)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/external/include)

link_directories(${PROJECT_BINARY_DIR}/external/lib)

add_subdirectory(src)


# add_executable(wavcav.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wavcav.cpp)
# target_link_libraries(wavcav.x cavity pedra wavcav pwl wem getkw ${LIBS})
# 
# add_executable(pcmsolver.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/pcmsolver.cpp)
# target_link_libraries(pcmsolver.x solver green cavity pedra metal pwl wem wavcav getkw ${LIBS})
# 
# add_executable(pwl.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/pwllib/main_pwl.c)
# target_link_libraries(pwl.x pwl wem ${LIBS})
# 
# add_executable(pwc.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wemlib/main.c)
# target_link_libraries(pwc.x wem ${LIBS})
# 
# add_executable(cav.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wavcavlib/main.c)
# target_link_libraries(cav.x wavcav ${LIBS})

