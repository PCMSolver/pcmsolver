cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (pcmsolver CXX C Fortran)

message("Fortran compiler in pcmsolver ${CMAKE_Fortran_COMPILER}")

set (pcmsolver_VERSION_MAJOR 0)
set (pcmsolver_VERSION_MINOR 1)
# build executable

option(PCMSOLVER_ENABLE_PCMGREEN "Use Green´s function library" ON)
option(PCMSOLVER_ENABLE_PCMCAV "Use Cavity generator library" ON)

include(ExternalProject)

set (ExternalProjectCMakeArgs  
  -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCBLAS_ROOT=${CBLAS_ROOT}
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


