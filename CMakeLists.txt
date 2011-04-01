cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (pcmsolver CXX C Fortran)

set (pcmsolver_VERSION_MAJOR 0)
set (pcmsolver_VERSION_MINOR 1)
# build executable
set (EXECUTABLE "pcmsolver.exe")

option(PCMSOLVER_ENABLE_PCMGREEN "Use Green´s function library" ON)

set (pcmsolver_SOURCES
pcmsolver.cpp
TraditionalPCMOperator.cpp 
TraditionalPCMOperator.h)

if (PCMSOLVER_ENABLE_XCFUN)
	ExternalProject_Add(pcmgreen
		PREFIX ${CMAKE_SOURCE_DIR}/external
		SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/pcmgreen
		BINARY_DIR ${CMAKE_SOURCE_DIR}/external/pcmgreen-build
		DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
		DOWNLOAD_COMMAND git submodule update
		INSTALL_DIR ${CMAKE_SOURCE_DIR}/external
		CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
		)
	set(HAVE_PCMGREEN TRUE)
endif()

ExternalProject_Add_Step(pcmgreen init
	COMMAND git submodule init
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

find_package(Eigen)

include_directories (${EIGEN_INCLUDE_DIRS})

#add_definitions(-I/home/notasoul/luca/eigen)

add_executable ("pcmsolver.exe" ${pcmsolver_SOURCES})

# install executables and scripts
install (TARGETS ${EXECUTABLE} 
         RUNTIME DESTINATION "bin")

install (PROGRAMS ${SCRIPTS}
         DESTINATION "bin") 



