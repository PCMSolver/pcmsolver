cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (pcmsolver CXX C Fortran)

cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set (pcmsolver_VERSION_MAJOR 0)
set (pcmsolver_VERSION_MINOR 1)
# build executable

option(PCMSOLVER_ENABLE_PCMGREEN "Use Green´s function library" ON)
option(PCMSOLVER_ENABLE_PCMCAV "Use Cavity generator library" ON)

include(SaveCompilerFlags)

If (NOT DEFINED HAVE_CXX_FLAGS)
if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-unused")
	if (ENABLE_CODE_COVERAGE)
		set (CMAKE_CXX_FLAGS 
			"${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
		set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
	endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
	set (CMAKE_CXX_FLAGS "-Wno-unknown-pragmas")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -debug -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-debug -O3 -DNDEBUG")
	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
endif ()
SaveCompilerFlags(CXX)
endif ()

include(ExternalProject)

set (ExternalProjectCMakeArgs  
  -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCBLAS_ROOT=${CBLAS_ROOT}
  -DEIGEN3_ROOT=${EIGEN3_ROOT}
  )

if (PCMSOLVER_ENABLE_PCMGREEN)
  ExternalProject_Add(pcmgreen
	PREFIX ${CMAKE_SOURCE_DIR}/external
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/pcm-green
	BINARY_DIR ${CMAKE_SOURCE_DIR}/external/pcm-green-build
	DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
	DOWNLOAD_COMMAND git submodule update
	INSTALL_DIR ${CMAKE_SOURCE_DIR}/external
	CMAKE_ARGS ${ExternalProjectCMakeArgs}
	)
  set(HAVE_PCMGREEN TRUE)
endif()

if (PCMSOLVER_ENABLE_PCMCAV)
	ExternalProject_Add(pcmcav
		DOWNLOAD_COMMAND git submodule update
		DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
		PREFIX ${CMAKE_SOURCE_DIR}/external
		SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/pcmcav
		BINARY_DIR ${CMAKE_SOURCE_DIR}/external/pcmcav-build
		INSTALL_DIR ${CMAKE_SOURCE_DIR}/external
		CMAKE_ARGS ${ExternalProjectCMakeArgs}
#		CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
		)
	set(HAVE_PCMCAV TRUE)
endif()

ExternalProject_Add_Step(pcmgreen init
	COMMAND git submodule init
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

ExternalProject_Add_Step(pcmcav init
	COMMAND git submodule init
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

#find_package(Eigen)
find_package(Eigen3)
find_package(CBLAS)

link_directories (${CMAKE_SOURCE_DIR}/lib)
link_directories (${CMAKE_SOURCE_DIR}/external/lib)

include_directories (${EIGEN3_INCLUDE_DIR})
#include_directories (${EIGEN_INCLUDE_DIRS})
include_directories (${CMAKE_SOURCE_DIR}/include)
include_directories (${CMAKE_SOURCE_DIR}/external/include)

add_subdirectory(src)


