cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (pcmsolver CXX C Fortran)

set (pcmsolver_VERSION_MAJOR 0)
set (pcmsolver_VERSION_MINOR 1)
# build executable
set (EXECUTABLE "pcmsolver.exe")

option(PCMSOLVER_ENABLE_PCMGREEN "Use Green´s function library" ON)
option(PCMSOLVER_ENABLE_PCMCAV "Use Cavity generator library" ON)

include(ExternalProject)

set (pcmsolver_SOURCES
pcmsolver.cpp
TraditionalPCMOperator.cpp 
TraditionalPCMOperator.h)

if (PCMSOLVER_ENABLE_PCMGREEN)
	ExternalProject_Add(pcmgreen
		PREFIX ${CMAKE_SOURCE_DIR}
		SOURCE_DIR ${CMAKE_SOURCE_DIR}/pcm-green
		BINARY_DIR ${CMAKE_SOURCE_DIR}/pcm-green-build
		DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
		DOWNLOAD_COMMAND git submodule update
		INSTALL_DIR ${CMAKE_SOURCE_DIR}
		CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}
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
		CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
		)
	set(HAVE_PCMCAV TRUE)
endif()

ExternalProject_Add_Step(pcmgreen init
	COMMAND git submodule init
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

ExternalProject_Add_Step(pcmcav init
	COMMAND git submodule init
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external
	)

find_package(Eigen)

link_directories (${CMAKE_SOURCE_DIR}/lib)

include_directories (${EIGEN_INCLUDE_DIRS})

include_directories (${CMAKE_SOURCE_DIR}/include)

#add_definitions(-I/home/notasoul/luca/eigen)

add_executable (pcmsolver.exe ${pcmsolver_SOURCES})

target_link_libraries (pcmsolver.exe green metalsphere) 

# install executables and scripts
install (TARGETS ${EXECUTABLE} 
         RUNTIME DESTINATION "bin")

install (PROGRAMS ${SCRIPTS}
         DESTINATION "bin") 



