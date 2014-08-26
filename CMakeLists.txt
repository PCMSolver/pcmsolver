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
option(ENABLE_MPI        "Enable MPI parallelization"                                OFF)
option(ENABLE_OPENMP     "Enable OpenMP parallelization"                             OFF)
option(ENABLE_TESTS      "Enable compilation of unit tests"                          OFF)
option(DISABLE_EIGEN_OWN "Do not use Eigen3 headers shipped with the module"         OFF)
option(ENABLE_EIGEN_MKL  "Enable Eigen3 automatic fallback to some of Intel MKL algorithms"                 OFF)

#
#  CMake modules
#
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/testing
    )

set(LIBS)
set(libs_to_merge)

include(ConfigVersion)
include(ConfigArchitecture)
include(CheckCXX11Features)
include(ConfigCompilerFlags)
include(ConfigDocumentation)
include(ConfigExternal)
include(ConfigSafeGuards)
include(GenericMacros)
include(ConfigGitRevision) # Has to come after ConfigVersion.cmake
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(include/FCMangle.hpp)
include(MergeStaticLibs)

find_package(PythonInterp REQUIRED)

set(Boost_USE_STATIC_LIBS   ON)
if(ENABLE_TESTS)
	find_package(Boost 1.46.0 COMPONENTS filesystem system unit_test_framework REQUIRED)
else()
	find_package(Boost 1.46.0 COMPONENTS filesystem system REQUIRED)
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
    ${PROJECT_BINARY_DIR}/modules
    )

set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
	-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
	-DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}
   )

add_external(libgetkw)

set(TaylorCMakeArgs "${ExternalProjectCMakeArgs}")
list(REMOVE_ITEM TaylorCMakeArgs "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}" "-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}" "-DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}")

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
        CMAKE_ARGS ${TaylorCMakeArgs}
	)
else()
ExternalProject_Add(libtaylor
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${TaylorProjectCMakeArgs}
	)
endif()

configure_file(
    ${PROJECT_SOURCE_DIR}/Config.hpp.in
    ${PROJECT_BINARY_DIR}/include/Config.hpp
    )

set(INST_LIBDIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/share/libgetkw)
configure_file(${PROJECT_SOURCE_DIR}/tools/pcmsolver.in pcmsolver)
file(COPY ${PROJECT_BINARY_DIR}/pcmsolver 
  DESTINATION ${PROJECT_BINARY_DIR}/bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

configure_file(${PROJECT_SOURCE_DIR}/tools/touch_cmakelists.py.in touch_cmakelists.py)
file(COPY ${PROJECT_BINARY_DIR}/touch_cmakelists.py 
  DESTINATION ${PROJECT_BINARY_DIR}/bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

configure_file(${PROJECT_SOURCE_DIR}/tools/update_gh-pages.py.in update_gh-pages.py)
file(COPY ${PROJECT_BINARY_DIR}/update_gh-pages.py 
  DESTINATION ${PROJECT_BINARY_DIR}/bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

include_directories(${EIGEN3_INCLUDE_DIR}
	            ${Boost_INCLUDE_DIRS}
		    ${Boost_QUADRATURE}
		    ${PROJECT_SOURCE_DIR}
		    ${PROJECT_SOURCE_DIR}/external
		    ${PROJECT_BINARY_DIR}/include
                    ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include
                    ${CMAKE_Fortran_MODULE_DIRECTORY})

add_subdirectory(src)
include_directories(${header_dir_list})

link_directories(${Boost_LIBRARY_DIRS}
                 ${EXTERNAL_PROJECT_INSTALL_PREFIX}/external/lib
	         ${EXTERNAL_PROJECT_INSTALL_PREFIX}/lib)
list(APPEND
     LIBS
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

        set(LIBS pthread z ${LIBS})

	setup_boosttest()

        enable_testing()
	add_subdirectory(tests) # This must come last!!
endif()

# How to add a debug executable. Shouldn't be necessary, we have unit tests...
#add_executable(debug_wavcav.x ${PROJECT_SOURCE_DIR}/src/debug_wavcav.cpp)
#target_link_libraries(debug_wavcav.x solver green cavity wavcav pwl wem utils ${LIBS})

# Merge static libs from subfolders into one static lib: libpcm.a
# This is the very last thing we do, i.e. DO NOT add anything depending on this target!!!
merge_static_libs(pcm
		  "${libs_to_merge}"
		 )
install(TARGETS pcm ARCHIVE DESTINATION lib)

