cmake_minimum_required (VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

message("${CMAKE_SOURCE_DIR}")

project (pcmsolver)
enable_language(CXX C Fortran)

set (pcmsolver_VERSION_MAJOR 0)
set (pcmsolver_VERSION_MINOR 1)
# build executable

# guard against in-source builds

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# Guard against bad build-type strings

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release (case-insensitive).")
endif()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/${CMAKE_BUILD_TYPE}")

include(SaveCompilerFlags)
include(FCompilers)
include(CCompilers)
include(CXXCompilers)

find_package(Git)
include(ExternalProject)

add_custom_command(
	OUTPUT ${PROJECT_SOURCE_DIR}/external/git-submodule-done
	COMMAND ${CMAKE_COMMAND} -E touch 
     	${PROJECT_SOURCE_DIR}/external/git-submodule-done
	COMMAND git submodule init
	COMMAND git submodule update
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	)

add_custom_target(git_submodules ALL
	DEPENDS ${PROJECT_SOURCE_DIR}/external/git-submodule-done
	)

add_custom_target(external)

set (ExternalProjectCMakeArgs  
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCBLAS_ROOT=${CBLAS_ROOT}
  -DEIGEN3_ROOT=${EIGEN3_ROOT}
  )

ExternalProject_Add(libtaylor
	DOWNLOAD_COMMAND svn checkout http://libtaylor.googlecode.com/svn/trunk/ ${PROJECT_SOURCE_DIR}/external/libtaylor
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor-build
	INSTALL_DIR ${PROJECT_SOURCE_DIR}/external
	CMAKE_ARGS ${ExternalProjectCMakeArgs}
	)

ExternalProject_Add(libgetkw
	DOWNLOAD_COMMAND git submodule update
	DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libgetkw
	BINARY_DIR ${PROJECT_SOURCE_DIR}/external/libgetkw-build
	INSTALL_DIR ${PROJECT_SOURCE_DIR}/external
	CMAKE_ARGS ${ExternalProjectCMakeArgs}
	)

add_dependencies(external libgetkw libtaylor)

include(Math)

find_package(Eigen3)

#cmake should exit gracefully if math libs are not found

add_definitions(-DSYS_LINUX)

#link_directories (${PROJECT_SOURCE_DIR}/lib)
link_directories (${PROJECT_SOURCE_DIR}/external/lib)

include_directories (${EIGEN3_INCLUDE_DIR})
include_directories (${PROJECT_SOURCE_DIR}/include)
include_directories (${PROJECT_SOURCE_DIR}/external/include)

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

