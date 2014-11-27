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

FortranCInterface_VERIFY(CXX)
init_FCMangle()
    
set(EMBEDDED_PYTHON FALSE)
# Try to configure Python, if not possible disable Python embedding 
if(ENABLE_PYTHON_EMBEDDING)
   # libutil and libdl are needed in case only the static version of libpython is available
   if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      find_package(Util REQUIRED)
      if(NOT LIBUTIL_FOUND)
         message(STATUS "Python embedding DISABLED (libutil NOT FOUND)")
      endif()	   
   endif()
   
   find_package(DL REQUIRED)
   if(NOT LIBDL_FOUND)
      message(STATUS "Python embedding DISABLED (libdl NOT FOUND)")
   endif()
   include(ConfigPython)
   # We set EMBEDDED_PYTHON to TRUE enabling the embedded Python code
   # only if the interpreter and the headers/libraries were found
   if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND) 
      set(EMBEDDED_PYTHON TRUE)
      add_definitions(-DEMBEDDED_PYTHON)
      include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
      message(STATUS "Python embedding ENABLED")
   else()
      message(STATUS "Python embedding DISABLED")
      message(STATUS "   Dependencies not satisfied")
   endif()
else(ENABLE_PYTHON_EMBEDDING)
   # The Python interpreter is anyway needed to run the MergeStaticLibs.py script
   find_package(PythonInterp 2.4 REQUIRED)
   message(STATUS "Python embedding DISABLED")
endif(ENABLE_PYTHON_EMBEDDING)

find_package(ZLIB REQUIRED) 
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
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
message(" : ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")

#  Eigen 3 stuff
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
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

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

# To be included as system headers
include_directories(SYSTEM 
		    ${PROJECT_SOURCE_DIR}/external
                    ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)
# To be included as project headers
include_directories(${PROJECT_SOURCE_DIR}
		    ${PROJECT_BINARY_DIR}/include
                    ${CMAKE_Fortran_MODULE_DIRECTORY})

add_subdirectory(src)
include_directories(${header_dir_list})

link_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/external/lib
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
