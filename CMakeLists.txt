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
option(ENABLE_TESTS            "Enable compilation of unit tests"                                 ON)
option(ENABLE_CODE_COVERAGE    "Enable code coverage"                                             OFF)
option(DISABLE_EIGEN_OWN       "Do not use Eigen3 headers shipped with the module"                OFF)
option(ENABLE_EIGEN_MKL        "Enable Eigen3 automatic fallback to some of Intel MKL algorithms" OFF)
option(ENABLE_64BIT_INTEGERS   "Enable 64-bit integers"                                           OFF)
option(ENABLE_CXX11_SUPPORT    "Enable C++11 compiler support"                                    ON)
# Developers only options
option(ENABLE_LOGGER           "Enable logger: needs C++11 support"                               ON)
option(ENABLE_TIMER            "Enable timer"                                                     ON)
#
#  CMake modules
#
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/build-info
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
include(ExternalProject)
include(GenericMacros)
include(MergeStaticLibs)
# Access to the fortran_enabler macro
include(FortranEnabler)
fortran_enabler()

# If C++11 support was not detected or was disabled by the user, disable the logger
get_property(has_cxx11_support GLOBAL PROPERTY CXX11_SUPPORT)
if(NOT has_cxx11_support)
   set(ENABLE_LOGGER OFF)
endif()

# The Python interpreter is needed to run the MergeStaticLibs.py script
# If the PYTHON_INTERPRETER variable is not empty, then the user passed a custom version
# of the interpreter to be used. If that's the case, set PYTHON_EXECUTABLE accordingly and
# proceed with detection of Python interpreter.
if("${PYTHON_INTERPRETER}" STREQUAL "")
   find_package(PythonInterp 2.4 REQUIRED)
else()
   set(PYTHONINTERP_FOUND TRUE)
   set(PYTHON_EXECUTABLE "${PYTHON_INTERPRETER}")
   find_package(PythonInterp 2.4 REQUIRED)
endif()

find_package(ZLIB REQUIRED)
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})
# librt   : needed for the timers (Linux only)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   find_package(RT REQUIRED)
   include_directories(SYSTEM "${LIBRT_INCLUDE_DIR}")
   link_libraries("${LIBRT_LIBRARIES}")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external)
   set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
else()
   set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
endif()

# Boost detection
include(ConfigBoost)
link_directories(${Boost_LIBRARY_DIRS})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

#  Eigen 3 stuff
if("${EIGEN3_ROOT}" STREQUAL "")
   if(DISABLE_EIGEN_OWN)
      find_package(Eigen3 3.0.0)
      message(STATUS "Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
   else()
      set(EIGEN3_ROOT ${CMAKE_SOURCE_DIR}/external/eigen3)
      set(EIGEN3_INCLUDE_DIR ${EIGEN3_ROOT}/include/eigen3)
      install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/eigen3 DESTINATION ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)
      message(STATUS "Eigen 3.2.0 is located here: " ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include/eigen3)
   endif()
else()
   find_package(Eigen3 3.0.0)
   message(STATUS "Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

if(ENABLE_EIGEN_MKL)
   message(STATUS "ENABLE_EIGEN_MKL option requires at least Eigen 3.1.0 and Intel MKL 10.3")
   message(STATUS "   Be sure you have read http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html")
   include(ConfigMath)
   set(EIGEN_USE_MKL_ALL ON)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_BINARY_DIR}/modules)

set(GetkwCMakeArgs
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
   -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
   -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}
   )

ExternalProject_Add(libgetkw
	PREFIX ${CMAKE_SOURCE_DIR}/external
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/libgetkw
	BINARY_DIR ${CMAKE_BINARY_DIR}/external/libgetkw-build
    STAMP_DIR ${CMAKE_BINARY_DIR}/external/libgetkw-stamp
    TMP_DIR ${CMAKE_BINARY_DIR}/external/libgetkw-tmp
    INSTALL_DIR ${CMAKE_BINARY_DIR}/external
    CMAKE_ARGS ${GetkwCMakeArgs})
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
	PREFIX ${CMAKE_SOURCE_DIR}/external
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${CMAKE_BINARY_DIR}/external/libtaylor-build
    STAMP_DIR ${CMAKE_BINARY_DIR}/external/libtaylor-stamp
    TMP_DIR ${CMAKE_BINARY_DIR}/external/libtaylor-tmp
    INSTALL_DIR ${CMAKE_BINARY_DIR}/external
    CMAKE_ARGS ${TaylorCMakeArgs})

configure_files()

# To be included as system headers
include_directories(SYSTEM
		    ${CMAKE_SOURCE_DIR}/external
                    ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)
# To be included as project headers
include_directories(${CMAKE_SOURCE_DIR}
		    ${CMAKE_BINARY_DIR}/include
                    ${CMAKE_Fortran_MODULE_DIRECTORY})

if(ENABLE_TESTS)
  # Build also standalone executables in src/bin
  set(BUILD_STANDALONE TRUE)
endif()
add_subdirectory(src)
include_directories(${header_dir_list})

link_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/external/lib
	         ${EXTERNAL_PROJECT_INSTALL_PREFIX}/lib)

install(FILES
	${CMAKE_BINARY_DIR}/bin/pcmsolver.py
	DESTINATION bin)

install(FILES
	${CMAKE_BINARY_DIR}/include/Config.hpp
	DESTINATION include)

install(FILES
	${CMAKE_BINARY_DIR}/include/FCMangle.hpp
	DESTINATION include)

install(FILES
	${CMAKE_BINARY_DIR}/include/Includer.hpp
	DESTINATION include)

if(ENABLE_TESTS)
   include(ConfigFramework)
   include(TestingMacros)

   setup_boosttest()

   enable_testing()
   add_subdirectory(tests) # This must come last!!
endif()

# How to add a debug executable. Shouldn't be necessary, we have unit tests...
#add_executable(debug_wavcav.x ${CMAKE_SOURCE_DIR}/src/debug_wavcav.cpp)
#target_link_libraries(debug_wavcav.x solver green cavity wavcav pwl wem utils ${LIBS})

# Merge static libs from subfolders into one static lib: libpcm.a
# This is the very last thing we do, i.e. DO NOT add anything depending on this target!!!
merge_static_libs(pcm "${libs_to_merge}")
install(TARGETS pcm ARCHIVE DESTINATION lib)

# configuration time in UTC
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import datetime; print(datetime.datetime.utcnow())"
    OUTPUT_VARIABLE CONFIGURATION_TIME
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

include(BuildInfo)
# Printout useful info
include(ConfigInfo)
