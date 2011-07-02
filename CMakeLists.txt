cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

message("${CMAKE_SOURCE_DIR}")

project (pcmsolver CXX C Fortran)

set (pcmsolver_VERSION_MAJOR 0)
set (pcmsolver_VERSION_MINOR 1)
# build executable

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/${CMAKE_BUILD_TYPE}")

include(Compilers)
include(SaveCompilerFlags)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

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

find_package(Git)
include(ExternalProject)

add_custom_command(
	OUTPUT ${CMAKE_SOURCE_DIR}/external/git-submodule-done
	COMMAND ${CMAKE_COMMAND} -E touch 
		${CMAKE_SOURCE_DIR}/external/git-submodule-done
	COMMAND git submodule init
	COMMAND git submodule update
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

add_custom_target(git_submodules ALL
	DEPENDS ${CMAKE_SOURCE_DIR}/external/git-submodule-done
	)

add_custom_target(external)

set (ExternalProjectCMakeArgs  
  -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCBLAS_ROOT=${CBLAS_ROOT}
  -DEIGEN3_ROOT=${EIGEN3_ROOT}
  )

ExternalProject_Add(libgetkw
	DOWNLOAD_COMMAND git submodule update
	DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
	PREFIX ${CMAKE_SOURCE_DIR}/external
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/libgetkw
	BINARY_DIR ${CMAKE_SOURCE_DIR}/external/libgetkw-build
	INSTALL_DIR ${CMAKE_SOURCE_DIR}/external
	CMAKE_ARGS ${ExternalProjectCMakeArgs}
	)

add_dependencies(external libgetkw)

find_package(Eigen3)
find_package(CBLAS)
find_package(BLAS)
find_package(LAPACK)
if(BLAS_FOUND AND LAPACK_FOUND)
    set(MATH_LIBS
        blas
        lapack
        )
    message("-- Set system default math libraries: ${MATH_LIBS} (${LAPACK_lapack_LIBRARY} ${BLAS_f77blas_LIBRARY})")
    set(MATH_FOUND TRUE)
endif()
set (LIBS ${MATH_LIBS})

add_definitions(-DSYS_LINUX)

link_directories (${CMAKE_SOURCE_DIR}/lib)
link_directories (${CMAKE_SOURCE_DIR}/external/lib)

include_directories (${EIGEN3_INCLUDE_DIR})
include_directories (${CMAKE_SOURCE_DIR}/include)
include_directories (${CMAKE_SOURCE_DIR}/external/include)

add_subdirectory(src)


