cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(PCMSolver)
enable_language(CXX C Fortran)

# There still is nothing going on inside the code regarding MPI, OpenMP and tests!!
option(ENABLE_MPI        "Enable MPI parallelization"                                OFF)
option(ENABLE_OPENMP     "Enable OpenMP parallelization"                             OFF)
option(ENABLE_TESTS      "Enable compilation of unit tests"                          OFF)
option(DISABLE_EIGEN_OWN "Do not use Eigen3 headers shipped with the module"         OFF)
option(ENABLE_EIGEN_MKL  "Enable Eigen3 automatic fallback to some of Intel MKL algorithms"                 OFF)

set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 0)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    )

set(LIBS)
set(SOURCES)

include(ConfigVersion)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(ConfigDocumentation)
include(ConfigExternal)
include(ConfigMath)
include(ConfigSafeGuards)
include(GenericMacros)
include(ConfigGitRevision)

find_package(PythonInterp REQUIRED)

find_package(Boost)

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
	        message("Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
	else()
		set(EIGEN3_ROOT ${PROJECT_SOURCE_DIR}/external/eigen3)
		set(EIGEN3_INCLUDE_DIR ${EIGEN3_ROOT}/include/eigen3)
		install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/eigen3 DESTINATION ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)
		message("Eigen 3.2.0 is located here: " ${EXTERNAL_PROJECT_INSTALL_PREFIX}/include/eigen3)
	endif()
else()
	find_package(Eigen3 3.0.0)
        message("Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
endif()

if(ENABLE_EIGEN_MKL)
	message("-- ENABLE_EIGEN_MKL option requires at least Eigen 3.1.0 and Intel MKL 10.3.\n Be sure you have read http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html")
        set(EIGEN_USE_MKL_ALL ON)
endif()


set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DPARENT_BINARY_DIR=${PROJECT_BINARY_DIR}
   )

add_external(libgetkw)

ExternalProject_Add(libtaylor
	SVN_REPOSITORY http://libtaylor.googlecode.com/svn/trunk/
       	DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
	PREFIX ${PROJECT_SOURCE_DIR}/external
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/libtaylor
	BINARY_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/libtaylor-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
	)

configure_file (
    ${PROJECT_SOURCE_DIR}/Config.hpp.in
    ${PROJECT_BINARY_DIR}/include/Config.hpp
    )

# ALL THIS STUFF NEEDS TO BE REVIEWED
set(INST_LIBDIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/share/libgetkw)
#set(PCMSOLVER_EXECUTABLE pcmsolver.x)
#set(PCMSOLVER_ERROR_LOG ${PROJECT_BINARY_DIR}/pcmsolver.err)
set(PCMSOLVER_EXECUTABLE ${PROJECT_BINARY_DIR}/pcmsolver.x)
option(PCMSOLVER_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
#ask jonas about path
configure_file(pcmsolver.in pcmsolver)
install(FILES ${PROJECT_BINARY_DIR}/pcmsolver 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ
	DESTINATION bin)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/include)

link_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/external/lib)
link_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/lib)

add_subdirectory(src)
include_directories(${header_dir_list})

add_library(pcm ${SOURCES})
install(TARGETS pcm ARCHIVE DESTINATION lib)

install(FILES 
	${PROJECT_BINARY_DIR}/include/Config.hpp 
	DESTINATION include)


if(ENABLE_TESTS)
	ExternalProject_Add(googletest
		SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk
	       	DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
		PREFIX ${PROJECT_SOURCE_DIR}/external
		SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/googletest
		BINARY_DIR ${PROJECT_BINARY_DIR}/external/googletest-build
	        STAMP_DIR ${PROJECT_BINARY_DIR}/external/googletest-stamp
        	TMP_DIR ${PROJECT_BINARY_DIR}/external/googletest-tmp
	        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        	CMAKE_ARGS ${ExternalProjectCMakeArgs}
		# Disable install step.
	        INSTALL_COMMAND ""
		)

        ExternalProject_Get_Property(googletest SOURCE_DIR)
	set(GTEST_INCLUDE_DIRS ${SOURCE_DIR}/include)
        ExternalProject_Get_Property(googletest BINARY_DIR)
	set(GTEST_LIBS_DIR ${BINARY_DIR})
	include_directories(${GTEST_INCLUDE_DIRS})

        enable_testing()
	file(GLOB_RECURSE TEST_SOURCE_FILES tests/*/*.cpp)
	foreach(child ${TEST_SOURCE_FILES})
		get_filename_component(test_name ${child} NAME_WE)
                add_executable(${test_name}.x ${child})      	
                add_dependencies(${test_name}.x googletest)
                target_link_libraries(${test_name}.x 
                interface solver green cavity pedra metal pwl wem wavcav utils getkw
                                      ${GTEST_LIBS_DIR}/libgtest.a
                                      ${GTEST_LIBS_DIR}/libgtest_main.a
	                              ${LIBS}
	                             )
                add_test(NAME ${test_name} COMMAND ${test_name}.x)
	endforeach()
endif()

# add_executable(wavcav.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wavcav.cpp)
# target_link_libraries(wavcav.x cavity pedra wavcav pwl wem getkw ${LIBS})
# 
# add_executable(pcmsolver.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/pcmsolver.cpp)
# target_link_libraries(pcmsolver.x solver green cavity pedra metal pwl wem wavcav getkw ${LIBS})
# 
# add_executable(pwl.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/pwllib/main_pwl.c)
# target_link_libraries(pwl.x pwl wem ${LIBS})
# 
# add_executable(pwc.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wemlib/main.c)
# target_link_libraries(pwc.x wem ${LIBS})
# 
# add_executable(cav.x EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/wavcavlib/main.c)
# target_link_libraries(cav.x wavcav ${LIBS})

