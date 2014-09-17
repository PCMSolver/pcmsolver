# Based on Ryan M. Richard work in PSI4
# First of all get BOOSTVER underscore separated
string(REGEX REPLACE "\\." "_" BOOSTVER ${BOOSTVER})
# Transform the ;-separated list to a ,-separated list (digested by the Boost build toolchain!)
string(REPLACE ";" "," needed_components "${needed_components}")

# Unpack Boost
add_custom_command(
    OUTPUT boost.unpacked 
    COMMAND tar -xzf ${PROJECT_SOURCE_DIR}/external/boost_${BOOSTVER}.tar.gz 
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external
    COMMENT "Unpacking Boost"
    DEPENDS ${PROJECT_SOURCE_DIR}/external/boost_${BOOSTVER}.tar.gz)

# We are going to build only what is necessary. We need to create a custom user-config.jam file
# To get boost to compile MPI we need to append "using mpi ;" to the end of the
# user-config.jam file. MPI_SENT will be the command we append
set(MPI_SENT "")
if(ENABLE_MPI)
	list(APPEND needed_components mpi)
        set(MPI_SENT "using mpi : ${CMAKE_CXX_COMPILER} \\;")
endif()
if(ENABLE_TESTS)
	# Replace unit_test_framework (used by CMake's find_package) with test (understood by Boost build toolchain)
	string(REPLACE "unit_test_framework" "test" needed_components "${needed_components}")
endif()
message(STATUS "  Libraries to be built: ${needed_components}")

# To get boost to be built with our compiler of choice we need these logic
# statements, then we write the options to user-config.jam
set(toolset "")
set(cxxflags "")
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  if(UNIX OR APPLE)
    set(toolset "intel-linux")
  else()
    set(toolset "intel-windows")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(toolset "clang")
else()
  if(APPLE)
    set(toolset "darwin")
  else()
   set(toolset "gcc")
  endif()
endif()
if(HAS_CXX11_SUPPORT)
	set(cxxflags "${cxxflags} ${CXX11_COMPILER_FLAGS}") 
endif()	  
message(STATUS "  Toolset to be used: ${toolset}")

# Write the user-config.jam file
file(WRITE  ${PROJECT_BINARY_DIR}/external/boost_${BOOSTVER}/user-config.jam "using ${toolset} : : ${CMAKE_CXX_COMPILER} \\;") 
file(APPEND ${PROJECT_BINARY_DIR}/external/boost_${BOOSTVER}/user-config.jam "${MPI_SENT}") 

# Run bootstrap.sh to configure the build. We will install in ${PROJECT_BINARY_DIR}/external/boost
add_custom_command(
    OUTPUT boost.configured
    COMMAND ./bootstrap.sh --with-toolset=${toolset}
            --with-libraries="${needed_components}"
	    --with-python=${PYTHON_EXECUTABLE}
	    --prefix=${PROJECT_BINARY_DIR}/external/boost 1> boost.configured.log 2> boost.configured.err
    COMMAND touch boost.configured	    
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external/boost_${BOOSTVER}
    DEPENDS boost.unpacked
    COMMENT "Configuring Boost")

# Append the MPI-related stuff to the project-config.jam file
file(APPEND ${PROJECT_BINARY_DIR}/external/boost_${BOOSTVER}/project-config.jam "${MPI_SENT}") 

# Build Boost
add_custom_command(
    OUTPUT boost.built
    COMMAND ./b2 1> boost.built.log 2> boost.built.err
    COMMAND touch boost.built
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external/boost_${BOOSTVER}
    DEPENDS boost.configured 
    COMMENT "Building Boost")

# Install Boost
add_custom_command(
    OUTPUT boost.installed
    COMMAND ./b2 install 1> boost.installed.log 2> boost.installed.err
    COMMAND touch boost.installed
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external/boost_${BOOSTVER}
    DEPENDS boost.built 
    COMMENT "Installing Boost")

add_custom_target(custom_boost ALL DEPENDS boost.installed)

# Set variables usually set by find_package by hand
set(Boost_FOUND TRUE)
set(Boost_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/external/boost/include) 
set(Boost_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/external/boost/lib)
# We will link statically, so just set the Boost_<C>_LIBRARY for the static library 
set(Boost_CHRONO_LIBRARY              ${PROJECT_BINARY_DIR}/external/boost/lib/libboost_chrono.a)
set(Boost_FILESYSTEM_LIBRARY          ${PROJECT_BINARY_DIR}/external/boost/lib/libboost_filesystem.a)
set(Boost_PYTHON_LIBRARY              ${PROJECT_BINARY_DIR}/external/boost/lib/libboost_python.a)
set(Boost_SYSTEM_LIBRARY              ${PROJECT_BINARY_DIR}/external/boost/lib/libboost_system.a)
set(Boost_TIMER_LIBRARY               ${PROJECT_BINARY_DIR}/external/boost/lib/libboost_timer.a)
if(ENABLE_TESTS)
	set(Boost_UNIT_TEST_FRAMEWORK_LIBRARY ${PROJECT_BINARY_DIR}/external/boost/lib/libboost_unit_test_framework.a)
endif()	
