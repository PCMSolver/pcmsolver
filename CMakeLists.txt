cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
# Set CMake policies here
if(CMAKE_VERSION VERSION_GREATER 2.8.3)
   # CMP0017 was introduced in CMake 2.8.4
   # Always prefer system modules to local ones
   cmake_policy(SET CMP0017 NEW)
endif()
#
#  Declare project name and programming languages
#
project(PCMSolver CXX C Fortran)
#  
#  Options  
#   There still is nothing going on inside the code regarding MPI and OpenMP!!
#
option(ENABLE_MPI              "Enable MPI parallelization"                                       OFF) 
option(ENABLE_OPENMP           "Enable OpenMP parallelization"                                    OFF)
option(ENABLE_VECTORIZATION    "Enable vectorization"                                             OFF)
option(ENABLE_TESTS            "Enable compilation of unit tests"                                 OFF)
option(ENABLE_CODE_COVERAGE    "Enable code coverage"                                             OFF)
option(DISABLE_EIGEN_OWN       "Do not use Eigen3 headers shipped with the module"                OFF)
option(ENABLE_EIGEN_MKL        "Enable Eigen3 automatic fallback to some of Intel MKL algorithms" OFF)
option(ENABLE_PYTHON_EMBEDDING "Enable Python embedding"                                          ON)
option(ENABLE_64BIT_INTEGERS   "Enable 64-bit integers"                                           OFF)
option(ENABLE_CXX11_SUPPORT    "Enable C++11 compiler support"                                    ON)
#
#  CMake modules
#
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/testing)

set(header_dir_list)
set(libs_to_merge)

include(ConfigSafeGuards)
include(ConfigVersion)
include(ConfigArchitecture)
include(CheckIncludeFile)
include(ConfigCompilerFlags)
include(PerTargetCompilerFlags)
include(ConfigDocumentation)
include(ConfigExternal)
include(GenericMacros)
include(ConfigGitRevision) # Has to come after ConfigVersion.cmake
include(MergeStaticLibs)
include(FortranCInterface)
include(CheckPython)

FortranCInterface_VERIFY(CXX)
init_FCMangle()
    
set(EMBEDDED_PYTHON FALSE)
# Check Python:
if(ENABLE_PYTHON_EMBEDDING)
   # 1. a Python interpreter of the suitable version exists
   if("${PYTHON_INTERPRETER}" STREQUAL "")
      find_package(PythonInterp 2.4 REQUIRED)
   else()
      set(PYTHONINTERP_FOUND TRUE)
      set(PYTHON_EXECUTABLE "${PYTHON_INTERPRETER}")
      get_Python_version_string(pyVersion)
      set(PYTHON_VERSION_STRING "${pyVersion}")
      message(STATUS "Passed PythonInterp: ${PYTHON_INTERPRETER} (found suitable version \"${PYTHON_VERSION_STRING}\", minimum required is \"2.4\")")
   endif()	
   # 2. Python libraries and headers of the same version exist
   # Set variables to help find Python library that is compatible with interpreter
   # Copied from https://bitbucket.org/fenics-project/dolfin
   if(PYTHONINTERP_FOUND)
      # Get Python include path from Python interpretter
      execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                            "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
                     OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
                     RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)
      # Get Python library path from interpreter
      execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                           "import os, sys, inspect; sys.stdout.write(os.path.split(os.path.split(inspect.getfile(inspect))[0])[0])"
                     OUTPUT_VARIABLE _PYTHON_LIB_PATH
                     RESULT_VARIABLE _PYTHON_LIB_RESULT)
      # Set include path, if returned by interpreter
      if("${_PYTHON_INCLUDE_RESULT}" STREQUAL "0")
         set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH})
      endif()
      # Add a search path for Python library based on output from
      # interpreter
      set(CMAKE_LIBRARY_PATH_SAVE ${CMAKE_LIBRARY_PATH})
      if("${_PYTHON_LIB_RESULT}" STREQUAL "0")
         set(CMAKE_LIBRARY_PATH ${_PYTHON_LIB_PATH})
      endif()
      # Find Pythons libs
      find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
      # Check headers
      include(CheckIncludeFiles)
      set(CMAKE_REQUIRED_INCLUDES "${PYTHON_INCLUDE_DIRS}")
      check_include_files(Python.h HAS_PYTHON_H)
      check_include_files(pyconfig.h HAS_PYCONFIG_H)
      if(NOT HAS_PYTHON_H)
	 message(STATUS "No Python.h header found!!")	   
      endif()	   
      if(NOT HAS_PYCONFIG_H)
	 message(STATUS "No pyconfig.h header found!!")	   
      endif()	   
      unset(CMAKE_REQUIRED_INCLUDES)
      # 3. that we can link against those libraries
      check_Python_compiles(pyCompiles)
      if(NOT pyCompiles)
	 message(STATUS "Cannot link against Python!!")
      endif()	 
   endif()
   # Iff the answer is "YES" to all of the above, we set
   # EMBEDDED_PYTHON to TRUE enabling the embedded Python code
   if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND 
      AND HAS_PYTHON_H AND HAS_PYCONFIG_H AND pyCompiles)
      set(EMBEDDED_PYTHON TRUE)
      add_definitions(-DEMBEDDED_PYTHON)
      message(STATUS "Python embedding ENABLED")
   else()
      message(STATUS "Python embedding DISABLED")
      message(STATUS "   Dependencies not satisfied")
   endif()
else()
   # The Python interpreter is anyway needed to run the MergeStaticLibs.py script
   find_package(PythonInterp 2.4 REQUIRED)
   message(STATUS "Python embedding DISABLED")
endif()
find_package(ZLIB REQUIRED) 

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
	set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
else()
	set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
endif()

# Just change the Boost version number here
set(BOOSTVER 1.54.0)
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
	set(CUSTOM_BOOST_LOCATION ${EXTERNAL_PROJECT_INSTALL_PREFIX}/boost)
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

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules)

set(GetkwCMakeArgs
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
   -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
   -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS})

add_external(libgetkw "" "${GetkwCMakeArgs}")
# In case we need to build Boost, let libgetkw depend on custom_boost
# so that's the first thing it's done
if(BUILD_CUSTOM_BOOST)
	add_dependencies(libgetkw custom_boost)
endif()
set(GETKW_PYTHON_DIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/share/libgetkw)

set(TaylorCMakeArgs 
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   )

ExternalProject_Add(libtaylor
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${TaylorCMakeArgs})

configure_files()

# To be includes as system headers
include_directories(SYSTEM 
	            ${EIGEN3_INCLUDE_DIR}
	            ${Boost_INCLUDE_DIRS}
		    ${PYTHON_INCLUDE_DIRS}
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
	         ${PYTHON_LIBRARIES}
                 ${EXTERNAL_PROJECT_INSTALL_PREFIX}/external/lib
	         ${EXTERNAL_PROJECT_INSTALL_PREFIX}/lib)

install(FILES 
	${PROJECT_BINARY_DIR}/bin/pcmsolver.py 
	DESTINATION bin)

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
	
	setup_boosttest()

        enable_testing()
	add_subdirectory(tests) # This must come last!!
endif()

# How to add a debug executable. Shouldn't be necessary, we have unit tests...
#add_executable(debug_wavcav.x ${PROJECT_SOURCE_DIR}/src/debug_wavcav.cpp)
#target_link_libraries(debug_wavcav.x solver green cavity wavcav pwl wem utils ${LIBS})

# Merge static libs from subfolders into one static lib: libpcm.a
# This is the very last thing we do, i.e. DO NOT add anything depending on this target!!!
merge_static_libs(pcm "${libs_to_merge}")
install(TARGETS pcm ARCHIVE DESTINATION lib)

# Printout useful info
include(ConfigInfo)
