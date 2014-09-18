cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#
#  Declare project name and programming languages
#
project(PCMSolver)
enable_language(CXX C Fortran)
#  
#  Options  
#   There still is nothing going on inside the code regarding MPI and OpenMP!!
#
option(ENABLE_MPI           "Enable MPI parallelization"                                       OFF)
option(ENABLE_OPENMP        "Enable OpenMP parallelization"                                    OFF)
option(ENABLE_VECTORIZATION "Enable vectorization"                                             OFF)
option(ENABLE_TESTS         "Enable compilation of unit tests"                                 OFF)
option(ENABLE_CODE_COVERAGE "Enable code coverage"                                             OFF)
option(DISABLE_EIGEN_OWN    "Do not use Eigen3 headers shipped with the module"                OFF)
option(ENABLE_EIGEN_MKL     "Enable Eigen3 automatic fallback to some of Intel MKL algorithms" OFF)
#
#  CMake modules
#
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/testing)

set(LIBS)
set(header_dir_list)
set(libs_to_merge)

include(ConfigVersion)
include(ConfigArchitecture)
include(CheckCXX11Features)
include(ConfigCompilerFlags)
include(PerTargetCompilerFlags)
include(ConfigDocumentation)
include(ConfigExternal)
include(ConfigSafeGuards)
include(GenericMacros)
include(ConfigGitRevision) # Has to come after ConfigVersion.cmake
include(MergeStaticLibs)
include(FortranCInterface)

FortranCInterface_VERIFY(CXX)
init_FCMangle()

find_package(PythonInterp REQUIRED)
find_package(ZLIB REQUIRED) 

# Just change the Boost version number here
set(BOOSTVER 1.55.0)
set(BUILD_CUSTOM_BOOST FALSE)
# List all components needed (except mpi and unit_test_framework) here.
# Components additionally required in PSI4: python, serialization, thread (Might be useful in the future?)
# mpi and unit_test_framework will be added afterwards, if needed.
list(APPEND needed_components filesystem system)
set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
if(ENABLE_TESTS)
	list(APPEND needed_components unit_test_framework)
	find_package(Boost ${BOOSTVER} COMPONENTS "${needed_components}")
else()
	find_package(Boost ${BOOSTVER} COMPONENTS "${needed_components}")
endif()
if(NOT Boost_FOUND)
        # Set also variables usually set by find_package
	message(STATUS "  Build custom Boost ${BOOSTVER}")
	set(BUILD_CUSTOM_BOOST TRUE)
	set(CUSTOM_BOOST_LOCATION ${PROJECT_BINARY_DIR}/external/boost)
	string(REGEX REPLACE "\\." "0" Boost_VERSION ${BOOSTVER})
        math(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
        math(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
        math(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
	set(Boost_LIB_VERSION ${Boost_MAJOR_VERSION}_${Boost_MINOR_VERSION})
	add_subdirectory(external)
	set(Boost_FOUND TRUE)
	set(Boost_LIBRARIES "")
	# Read documentation in FindBoost.cmake for the difference between the singular and plural forms
	set(Boost_INCLUDE_DIR  ${CUSTOM_BOOST_LOCATION}/include)
	set(Boost_INCLUDE_DIRS ${CUSTOM_BOOST_LOCATION}/include) 
	set(Boost_LIBRARY_DIR  ${CUSTOM_BOOST_LOCATION}/lib)
	set(Boost_LIBRARY_DIRS ${CUSTOM_BOOST_LOCATION}/lib)
        # We will link statically, so just set the Boost_<C>_LIBRARY for the static library 
	foreach(_component ${needed_components})
		string(TOUPPER ${_component} _COMP)
		set(Boost_${_COMP}_FOUND TRUE)
		set(Boost_${_COMP}_LIBRARY libboost_${_component}-${Boost_LIB_VERSION}.a)
		list(APPEND Boost_LIBRARIES Boost_${_COMP}_LIBRARY)
	endforeach()
	# Set environment variables to avoid problems with libgetkw not recognizing Boost paths
	set(ENV{BOOST_INCLUDEDIR} ${Boost_INCLUDE_DIRS})
	set(ENV{BOOST_LIBRARYDIR} ${Boost_LIBRARY_DIRS})
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
	set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
else()
	set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
endif()

#
#  Eigen 3 stuff
#
if("${EIGEN3_ROOT}" STREQUAL "")
	if(DISABLE_EIGEN_OWN)
		find_package(Eigen3 3.0.0)
	        message(STATUS "Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
	else()
		set(EIGEN3_ROOT ${PROJECT_SOURCE_DIR}/external/eigen3)
		set(EIGEN3_INCLUDE_DIR ${EIGEN3_ROOT}/include/eigen3)
		install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/eigen3 DESTINATION ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)
		message(STATUS "Eigen 3.2.0 is located here: " ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include/eigen3)
	endif()
else()
	find_package(Eigen3 3.0.0)
        message(STATUS "Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
endif()

if(ENABLE_EIGEN_MKL)
	message(STATUS "ENABLE_EIGEN_MKL option requires at least Eigen 3.1.0 and Intel MKL 10.3") 
        message(STATUS "   Be sure you have read http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html")
        include(ConfigMath)
        set(EIGEN_USE_MKL_ALL ON)
endif()

#
#  Boost Numeric Quadrature stuff
#
set(Boost_QUADRATURE ${PROJECT_SOURCE_DIR}/external/boost)
install(DIRECTORY ${Boost_QUADRATURE} DESTINATION ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules)

set(ExternalProjectCMakeArgs
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
   -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
   -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS})

add_external(libgetkw)
# In case we need to build Boost, let libgetkw depend on custom_boost
# so that's the first thing it's done
if(BUILD_CUSTOM_BOOST)
	add_dependencies(libgetkw custom_boost)
endif()
set(INST_LIBDIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/share/libgetkw)

set(TaylorCMakeArgs "${ExternalProjectCMakeArgs}")
list(REMOVE_ITEM TaylorCMakeArgs 
     "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}" "-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}" "-DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}")
list(APPEND TaylorCMakeArgs -Wno-dev)

if(DEVELOPMENT_CODE)
ExternalProject_Add(libtaylor
	SVN_REPOSITORY http://libtaylor.googlecode.com/svn/trunk/
       	DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${TaylorCMakeArgs})
else()
ExternalProject_Add(libtaylor
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${TaylorCMakeArgs})
endif()

configure_files()

# To be includes as system headers
include_directories(SYSTEM 
	            ${EIGEN3_INCLUDE_DIR}
	            ${Boost_INCLUDE_DIRS}
		    ${Boost_QUADRATURE}
		    ${ZLIB_INCLUDE_DIRS}
		    ${PROJECT_SOURCE_DIR}/external
                    ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)
# To be included as project headers
include_directories(${PROJECT_SOURCE_DIR}
		    ${PROJECT_BINARY_DIR}/include
                    ${CMAKE_Fortran_MODULE_DIRECTORY})

add_subdirectory(src)
include_directories(${header_dir_list})

link_directories(${Boost_LIBRARY_DIRS}
                 ${EXTERNAL_PROJECT_INSTALL_PREFIX}/external/lib
	         ${EXTERNAL_PROJECT_INSTALL_PREFIX}/lib)
list(APPEND LIBS
    ${Boost_FILESYSTEM_LIBRARY} 
    ${Boost_SYSTEM_LIBRARY})

install(FILES 
	${PROJECT_BINARY_DIR}/include/Config.hpp 
	DESTINATION include)

install(FILES 
	${PROJECT_BINARY_DIR}/include/FCMangle.hpp 
	DESTINATION include)

install(FILES 
	${PROJECT_BINARY_DIR}/include/Includer.hpp 
	DESTINATION include)

if(ENABLE_TESTS)
	include(ConfigFramework)
	include(TestingMacros)
	
	find_package(Threads REQUIRED)

	list(APPEND LIBS ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES})

	setup_boosttest()

        enable_testing()
	add_subdirectory(tests) # This must come last!!
endif()

# Merge static libs from subfolders into one static lib: libpcm.a
# This is the very last thing we do, i.e. DO NOT add anything depending on this target!!!
merge_static_libs(pcm "${libs_to_merge}")
install(TARGETS pcm ARCHIVE DESTINATION lib)
