if(NOT DEFINED ENV{CXXFLAGS})
  if(CMAKE_CXX_COMPILER_ID MATCHES Intel)

    # We require C++11 support from the compiler and standard library.
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.1")
      message(FATAL_ERROR "ICPC version must be at least 2015.0.1!")
    endif()

    set(_testfl ${CMAKE_BINARY_DIR}/test_gcc_version.cc)
    file(WRITE  ${_testfl} "
    #include <stdio.h>

    int main() {
        #ifdef __clang__
        printf(\"%d.%d.%d\", __clang_major__, __clang_minor__, __clang_patchlevel__);
        #else
        printf(\"%d.%d.%d\", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
        #endif
        return 0;
    }
    ")
    try_run(GCCV_COMPILES
            GCCV_RUNS
            ${CMAKE_BINARY_DIR} ${_testfl}
            RUN_OUTPUT_VARIABLE GCC_VERSION)
    message(STATUS "Found base compiler version ${GCC_VERSION}")
    file(REMOVE ${_testfl})

    if (APPLE)
        if (${GCC_VERSION} VERSION_LESS 6.1)
            message(FATAL_ERROR "${BoldYellow}Intel ICPC makes use of Clang (detected: ${GCC_VERSION}; required for C++11: 6.1) so this build won't work without CLANG intervention.\n${ColourReset}")
        endif()
    else ()
        if (${GCC_VERSION} VERSION_LESS 4.8)
            message(FATAL_ERROR "${BoldYellow}Intel ICPC makes use of GCC (detected: ${GCC_VERSION}; required for C++11: 4.9) so this build won't work without GCC intervention.\n${ColourReset}")
        endif()
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -debug -DDEBUG -Wall -Wuninitialized -Wno-unknown-pragmas -Wno-sign-compare")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  endif()
endif()
